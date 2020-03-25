#[macro_use]
extern crate log;
extern crate bio_utils;
extern crate edlib_sys;
extern crate env_logger;
extern crate last_tiling;
extern crate md5;
extern crate nalgebra as na;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
extern crate serde_json;
#[macro_use]
extern crate serde;
pub use find_breakpoint::critical_regions;
use rayon::prelude::*;
pub mod bipartite_matching;
mod eread;
pub mod find_breakpoint;
mod find_union;
pub mod utils;
pub use eread::*;
pub use find_breakpoint::initial_clusters;
pub use find_breakpoint::Cluster;
use find_breakpoint::ReadClassify;
use last_tiling::repeat::RepeatPairs;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::{HashMap, HashSet};
pub mod d3_data;
pub mod error_profile;
pub mod poa_clustering;
pub use poa_clustering::soft_clustering_poa;
use poa_hmm::Config;
pub mod variant_calling;
mod digamma;
// 200 * 100 = 20_000
const WINDOW_SIZE: usize = 200;
const OVERLAP: usize = 50;
const MIN_LEN: usize = 5_000;
const CONNECTION_THR: f64 = 20.;
type Read = Vec<(usize, Vec<u8>)>;
/// Main method. Decomposing the reads.
/// You should call "merge" method separatly(?) -- should be integrated with this function.
// TODO: Make a procedure to remove chimeric reads.
// -> If a read has 'junction', and in the "unassigned" bucket,
// it means such read can be safely removed.
pub fn decompose(
    encoded_reads: Vec<ERead>,
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
    config: &Config,
    answer: &HashMap<String, u8>,
    cluster_num: usize,
) -> Vec<(String, u8)> {
    let datasize = encoded_reads.len();
    let mut unassigned_reads: Vec<_> = vec![];
    let mut assigned_reads: Vec<_> = vec![];
    let mut labels: Vec<_> = vec![];
    for read in encoded_reads {
        let matched_cluster = initial_clusters
            .iter()
            .filter_map(|cr| if cr.has(read.id()) { Some(cr.id) } else { None })
            .nth(0);
        if let Some(idx) = matched_cluster {
            assigned_reads.push(read);
            labels.push(idx as u8);
        } else {
            unassigned_reads.push(read);
        }
    }
    assert_eq!(labels.len(), assigned_reads.len());
    assert_eq!(assigned_reads.len() + unassigned_reads.len(), datasize);
    let masked_region = get_masked_region(&initial_clusters, &contigs, repeats);
    // Remove contained reads.
    let mut forbidden = vec![];
    let assigned_reads: Vec<_> = assigned_reads
        .into_iter()
        .map(|mut read| {
            forbidden.push(vec![]);
            let seq = read
                .seq()
                .iter()
                .filter(|unit| !masked_region[unit.contig as usize][unit.unit as usize])
                .cloned()
                .collect();
            *read.seq_mut() = seq;
            read
        })
        .filter(|read| !read.is_empty())
        .collect();
    // This should not decrease the number of assigned reads,
    // as each assigned reads should escaped from repetitive regions.
    assert_eq!(assigned_reads.len(), labels.len());
    debug!(
        "Unassigned reads before removing contained:{}",
        unassigned_reads.len()
    );
    let unassigned_reads: Vec<_> = unassigned_reads
        .into_iter()
        .filter_map(|mut read| {
            // let crs: Vec<_> = initial_clusters
            //     .iter()
            //     .filter(|c| c.is_spanned_by(&read))
            //     .map(|c| c.id as u8)
            //     .collect();
            let seq: Vec<_> = read
                .seq()
                .iter()
                .filter(|unit| !masked_region[unit.contig()][unit.unit()])
                .cloned()
                .collect();
            if seq.is_empty() {
                debug!("Read {} would be an empty read.", read);
            }
            *read.seq_mut() = seq;
            if !read.is_empty() {
                // forbidden.push(crs);
                Some(read)
            } else {
                None
            }
        })
        .collect();
    debug!(
        "Unassigned reads after removing contained:{}",
        unassigned_reads.len()
    );
    let answer: Vec<_> = unassigned_reads
        .iter()
        .map(|read| match answer.get(read.id()) {
            Some(cl) => *cl,
            None => 0,
        })
        .collect();
    let dataset = vec![assigned_reads, unassigned_reads].concat();
    assert_eq!(forbidden.len(), dataset.len());
    let total_units = dataset.iter().map(|read| read.seq().len()).sum::<usize>();
    debug!(
        "There are {} reads and {} units.",
        dataset.len(),
        total_units
    );
    assert_eq!(answer.len() + labels.len(), dataset.len());
    let contigs: Vec<_> = (0..contigs.get_num_of_contigs())
        .map(|e| contigs.get_last_unit(e as u16).unwrap() as usize + 1)
        .collect();
    let predicts = clustering_chunking(
        &dataset,
        &labels,
        &forbidden,
        initial_clusters,
        cluster_num,
        &contigs,
        &answer,
        &config,
    );
    dataset
        .into_iter()
        .zip(predicts.iter())
        .map(|(read, cl)| (read.id().to_string(), *cl))
        .collect()
}

fn get_masked_region(
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
) -> Vec<Vec<bool>> {
    let mut masked: Vec<Vec<_>> = contigs
        .get_last_units()
        .into_iter()
        .map(|len| vec![false; len as usize + 1])
        .collect();
    let repeats = repeats.iter().flat_map(|repeat| {
        repeat.inner().iter().map(|rep| {
            let contig = rep.id();
            let s = rep.start_in_unit() as i32;
            let t = rep.end_in_unit() as i32;
            (contig, (s, t))
        })
    });
    let ranges: Vec<_> = initial_clusters
        .iter()
        .flat_map(|e| e.ranges())
        .chain(repeats)
        .collect();
    for (c, (s, t)) in ranges {
        debug!("Masking {}:{}-{}", c, s, t);
        let (s, t) = (s as usize, t as usize);
        let c = c as usize;
        if s <= t {
            masked[c][s..t].iter_mut().for_each(|e| *e = true);
        } else {
            masked[c][s..].iter_mut().for_each(|e| *e = true);
            masked[c][..t].iter_mut().for_each(|e| *e = true);
        }
    }
    masked
}

pub fn clustering_via_alignment(
    reads: &[usize],
    label: &[u8],
    forbidden: &[Vec<u8>],
    similarity: &[HashMap<usize, i64>],
    cluster_num: usize,
) -> Vec<u8> {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(reads.len() as u64);
    let mut predictions: Vec<_> = reads
        .iter()
        .map(|_| rng.gen_range(0, cluster_num))
        .map(|e| e as u8)
        .collect();
    for (idx, &l) in label.iter().enumerate() {
        predictions[idx] = l;
    }
    let allowed: Vec<_> = reads
        .iter()
        .map(|&read| {
            let mut al = vec![true; cluster_num];
            for &f in &forbidden[read] {
                al[f as usize] = false;
            }
            al
        })
        .collect();
    let mut clusters: Vec<Vec<usize>> = (0..cluster_num)
        .map(|cl| {
            predictions
                .iter()
                .enumerate()
                .filter_map(|(idx, &assign)| if assign == cl as u8 { Some(idx) } else { None })
                .collect::<Vec<_>>()
        })
        .collect();
    let mut is_updated = true;
    let border = label.len();
    while is_updated {
        is_updated = predictions
            .par_iter_mut()
            .zip(reads.par_iter())
            .skip(border)
            .map(|(pred, &target)| {
                let (assign, _) = clusters
                    .iter()
                    .zip(allowed[target].iter())
                    .enumerate()
                    .filter_map(
                        |(i, (cl, &is_allowed))| {
                            if is_allowed {
                                Some((i, cl))
                            } else {
                                None
                            }
                        },
                    )
                    .map(|(idx, cluster)| {
                        let mean_sim = cluster
                            .iter()
                            .map(|query| similarity[target][query])
                            .sum::<i64>()
                            / cluster.len() as i64;
                        (idx, mean_sim)
                    })
                    .max_by_key(|e| e.1)
                    .unwrap_or((0, -1));
                let assign = assign as u8;
                let is_updated = assign != *pred;
                *pred = assign;
                is_updated
            })
            .reduce(|| false, |p, q| p | q);
        clusters = (0..cluster_num)
            .map(|cl| {
                predictions
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, &assign)| if assign == cl as u8 { Some(idx) } else { None })
                    .collect::<Vec<_>>()
            })
            .collect();
    }
    predictions
}

fn create_windows(idx: usize, len: usize) -> Vec<(u16, u16, u16)> {
    (0..)
        .map(|i| i * (WINDOW_SIZE - OVERLAP))
        .take_while(|&start_pos| start_pos + WINDOW_SIZE / 2 < len || start_pos == 0)
        .map(|start_pos| {
            if start_pos + (WINDOW_SIZE - OVERLAP) + WINDOW_SIZE / 2 < len {
                (idx, start_pos, start_pos + WINDOW_SIZE)
            } else {
                (idx, start_pos, len)
            }
        })
        .map(|(x, y, z)| (x as u16, y as u16, z as u16))
        .collect()
}

fn select_within(
    (contig, start, end): (u16, u16, u16),
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    answer: &[u8],
) -> (Vec<ERead>, Vec<u8>, Vec<Vec<u8>>, Vec<u8>) {
    assert_eq!(data.len(), label.len() + answer.len());
    let (mut s_data, mut s_label, mut s_forbid, mut s_answer) = (vec![], vec![], vec![], vec![]);
    debug!("Selecting {}\t{}\t{}...", contig, start, end);
    let border = label.len();
    for i in 0..data.len() {
        let read = &data[i];
        let count = read.include_units(contig, start, end);
        if count > MIN_LEN / last_tiling::UNIT_SIZE {
            s_data.push(read.clone_within(contig, start, end));
            s_forbid.push(forbidden[i].clone());
            if i < border {
                s_label.push(label[i]);
            } else {
                s_answer.push(answer[i - border]);
            }
        }
    }
    (s_data, s_label, s_forbid, s_answer)
}

fn find_matching(prev: &[HashSet<String>], after: &[HashSet<String>]) -> Vec<(usize, usize)> {
    let graph: Vec<Vec<(usize, f64)>> = prev
        .iter()
        .map(|cl1| {
            after
                .iter()
                .enumerate()
                .map(|(idx, cl2)| (idx, sim(cl1, cl2)))
                .collect()
        })
        .collect();
    debug!("Dump graph");
    for (idx, edges) in graph.iter().enumerate() {
        for (to, weight) in edges.iter() {
            debug!("{}->({:.3})->{}", idx, weight, to);
        }
    }
    let res: Vec<(usize, usize)> = graph
        .into_iter()
        .enumerate()
        .flat_map(|(cl1, cl1_edges)| {
            cl1_edges
                .into_iter()
                .filter(|&(_, sim)| sim > CONNECTION_THR)
                .map(|(cl2, _)| (cl1, cl2))
                .collect::<Vec<(usize, usize)>>()
        })
        .collect();
    debug!("Path Selected.");
    for &(from, to) in &res {
        debug!("{}-({})-{}", from, sim(&prev[from], &after[to]), to);
    }
    res
}

// Define similarity between two cluster.
fn sim(a: &HashSet<String>, b: &HashSet<String>) -> f64 {
    a.intersection(&b).count() as f64
}

/// Clustering after chunking the reference into several chunks.
pub fn clustering_chunking(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    initial_clusters: &[Cluster],
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
    c: &Config,
) -> Vec<u8> {
    let windows: Vec<_> = contigs
        .iter()
        .enumerate()
        .flat_map(|(idx, len)| create_windows(idx, *len))
        .collect();
    debug!("The contig lengths:{:?}", contigs);
    debug!("There are {} windows in total.", windows.len());
    let clusterings: Vec<Vec<HashSet<String>>> = windows
        .iter()
        .map(|&region| {
            {
                let (contig, start, end) = region;
                debug!("{}:{}-{}", contig, start, end);
            }
            // Determine the number of the cluster.
            let cluster_num = initial_clusters
                .iter()
                .filter(|cl| cl.overlap(region))
                .count()
                .max(cluster_num - 1)
                + 1;
            debug!("Number of clusters:{}", cluster_num);
            let (data, label, forbidden, answer) =
                select_within(region, data, label, forbidden, answer);
            debug!("Number of Reads:{}", data.len());
            assert_eq!(data.len(), forbidden.len());
            assert_eq!(data.len(), label.len() + answer.len());
            let predictions = {
                let (da, la, fo, an) = (&data, &label, &forbidden, &answer);
                clustering(da, la, fo, cluster_num, an, c)
            };
            assert_eq!(predictions.len(), data.len());
            (0..cluster_num)
                .map(|cluster_idx| {
                    let cluster_idx = cluster_idx as u8;
                    predictions
                        .iter()
                        .zip(data.iter())
                        .filter_map(|(&p, r)| if p == cluster_idx { Some(r.id()) } else { None })
                        .map(|id| id.to_string())
                        .collect()
                })
                .collect()
        })
        .collect();
    // Merging.
    let mut fu = find_union::FindUnion::new(cluster_num * windows.len());
    for idx in 0..clusterings.len() - 1 {
        let prev_idx = idx;
        let after_idx = (idx + 1) % clusterings.len();
        let prev = &clusterings[prev_idx];
        let after = &clusterings[after_idx];
        assert!(prev.len() == cluster_num);
        assert!(after.len() == cluster_num);
        for (i, j) in find_matching(prev, after) {
            let i = i + prev_idx * cluster_num;
            let j = j + after_idx * cluster_num;
            fu.unite(i, j).unwrap();
        }
    }
    // Then, iteratively take components.
    let mut result: HashMap<String, u8> = HashMap::new();
    let mut current = 0;
    for i in 0..(cluster_num * windows.len()) {
        let parent = fu.find(i).unwrap();
        if parent != i {
            continue;
        }
        info!("Find cluster");
        for (w, clusters) in clusterings.iter().enumerate() {
            for (j, cluster) in clusters.iter().enumerate() {
                if parent == fu.find(j + w * cluster_num).unwrap() {
                    info!("{}:{}", w, j);
                    for id in cluster {
                        result.insert(id.clone(), current);
                    }
                }
            }
        }
        debug!("These clusters assigned as {}", current);
        current += 1;
    }
    data.iter()
        .map(|read| match result.get(read.id()) {
            Some(res) => *res,
            None => {
                info!("Read {} does not belong to any cluster.", read.id());
                0
            }
        })
        .collect()
}

pub fn clustering(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    answer: &[u8],
    c: &Config,
) -> Vec<u8> {
    assert_eq!(forbidden.len(), data.len());
    assert_eq!(label.len() + answer.len(), data.len());
    use poa_clustering::DEFAULT_ALN;
    let weights = soft_clustering_poa(data, label, forbidden, cluster_num, answer, c, &DEFAULT_ALN);
    //     variational_bayes_poa(data, label, forbidden, cluster_num, answer, c, &DEFAULT_ALN);
    debug!("WEIGHTS\tPrediction. Dump weights");
    assert_eq!(weights.len(), label.len() + answer.len());
    for (weight, ans) in weights.iter().zip(label.iter().chain(answer.iter())) {
        let weights: String = weight
            .iter()
            .map(|e| format!("{:.1}\t", e))
            .fold(String::new(), |x, y| x + &y);
        debug!("WEIGHTS\t{}{}", weights, ans);
    }
    weights
        .iter()
        .map(|weight| {
            assert_eq!(weight.len(), cluster_num);
            let (cl, _max): (u8, f64) = weight.iter().enumerate().fold(
                (0, -1.),
                |(i, m), (j, &w)| if m < w { (j as u8, w) } else { (i, m) },
            );
            cl
        })
        .collect()
}
