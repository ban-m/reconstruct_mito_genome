#[macro_use]
extern crate log;
extern crate bio_utils;
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
pub mod annotate_contigs_to_reference;
pub mod d3_data;
pub mod error_profile;
pub mod poa_clustering;
pub use poa_clustering::soft_clustering_poa;
use poa_hmm::Config;
mod digamma;
pub mod variant_calling;
// 200 * 100 = 25_000
const WINDOW_SIZE: usize = 250;
const OVERLAP: usize = 50;
const MIN_LEN: usize = 5_000;
const CONNECTION_THR: f64 = 10.;
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
    limit: u64,
) -> Vec<(String, u8)> {
    let datasize = encoded_reads.len();
    let masked_region = get_masked_region(&initial_clusters, &contigs, repeats);
    let mut forbidden = vec![];
    let mut labels: Vec<_> = vec![];
    let (assigned_reads, unassigned_reads): (Vec<_>, Vec<_>) = encoded_reads
        .into_iter()
        .map(|mut read| {
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
        .partition(|read| {
            let matched_cluster = initial_clusters
                .iter()
                .filter_map(|cr| if cr.has(read.id()) { Some(cr.id) } else { None })
                .nth(0);
            forbidden.push(vec![]);
            if let Some(idx) = matched_cluster {
                labels.push(idx as u8);
                true
            } else {
                false
            }
        });
    assert_eq!(labels.len(), assigned_reads.len());
    assert!(assigned_reads.len() + unassigned_reads.len() <= datasize);
    debug!(
        "Unassigned reads before removing contained:{}",
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
    // assert_eq!(forbidden.len(), dataset.len());
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
    let windows: Vec<_> = contigs
        .iter()
        .zip(masked_region.iter())
        .enumerate()
        .flat_map(|(idx, (len, mask))| create_windows(idx, *len, mask))
        .collect();
    let predicts = clustering_chunking(
        &dataset,
        (&labels, &answer),
        &forbidden,
        initial_clusters,
        cluster_num,
        &contigs,
        &config,
        limit,
        &windows,
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

fn create_windows(idx: usize, len: usize, mask: &[bool]) -> Vec<(u16, u16, u16)> {
    let mut windows = vec![];
    let mut start = 0;
    while start < len {
        // determine end position.
        let (mut end, mut unmasked_count) = (start, 0);
        while unmasked_count < WINDOW_SIZE && end < len {
            unmasked_count += if !mask[end] { 1 } else { 0 };
            end += 1;
        }
        if end + WINDOW_SIZE / 2 - OVERLAP < len {
            windows.push((idx as u16, start as u16, end as u16));
            start = end;
        } else {
            windows.push((idx as u16, start as u16, len as u16));
            break;
        }
    }
    windows
}

fn select_within(
    (contig, start, end): (u16, u16, u16),
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    answer: &[u8],
) -> ((Vec<ERead>, Vec<ERead>), Vec<u8>, Vec<Vec<u8>>, Vec<u8>) {
    assert_eq!(data.len(), label.len() + answer.len());
    let (mut s_data, mut s_label, mut s_forbid, mut s_answer) = (vec![], vec![], vec![], vec![]);
    debug!("Selecting {}\t{}\t{}...", contig, start, end);
    let mut shorter_reads = vec![];
    let border = label.len();
    for (idx, read) in data.iter().enumerate() {
        let original_len = read.seq().len();
        let read = read.clone_within(contig, start, end);
        if read.seq().len() > MIN_LEN / last_tiling::UNIT_SIZE {
            s_data.push(read);
            s_forbid.push(forbidden[idx].clone());
            if idx < border {
                s_label.push(label[idx]);
            } else {
                s_answer.push(answer[idx - border]);
            }
        } else if read.seq().len() > original_len / 2 {
            shorter_reads.push(read);
        }
    }
    ((s_data, shorter_reads), s_label, s_forbid, s_answer)
}

fn find_matching(prev: &[HashSet<String>], after: &[HashSet<String>]) -> Vec<(usize, usize)> {
    let node1 = prev.len();
    let node2 = after.len();
    debug!("Dump graph");
    let graph: Vec<Vec<(usize, f64)>> = prev
        .iter()
        .enumerate()
        .map(|(from, cl1)| {
            after
                .iter()
                .enumerate()
                .filter_map(|(to, cl2)| {
                    let sim = sim(cl1, cl2);
                    debug!("{}->({:.3})->{}", from, sim, to);
                    if sim > CONNECTION_THR {
                        Some((to, sim))
                    } else {
                        None
                    }
                })
                .collect()
        })
        .collect();
    let res = bipartite_matching::maximum_weight_matching(node1, node2, &graph);
    // let res: Vec<(usize, usize)> = graph
    //     .into_iter()
    //     .enumerate()
    //     .flat_map(|(cl1, cl1_edges)| {
    //         cl1_edges
    //             .into_iter()
    //             .filter(|&(_, sim)| sim > CONNECTION_THR)
    //             .map(|(cl2, _)| (cl1, cl2))
    //             .collect::<Vec<(usize, usize)>>()
    //     })
    //     .collect();
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
    (label, answer): (&[u8], &[u8]),
    forbidden: &[Vec<u8>],
    initial_clusters: &[Cluster],
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
    limit: u64,
    windows: &[(u16, u16, u16)],
) -> Vec<u8> {
    debug!("The contig lengths:{:?}", contigs);
    debug!("There are {} windows in total.", windows.len());
    for w in windows {
        debug!("{:?}", w);
    }
    let clusterings: Vec<Vec<HashSet<String>>> = windows
        .iter()
        .enumerate()
        .map(|(idx, &region)| {
            {
                let (contig, start, end) = region;
                debug!("{}/{}:{}-{}", idx, contig, start, end);
            }
            // Determine the number of the cluster.
            let cluster_num = initial_clusters
                .iter()
                .filter(|cl| cl.overlap(region))
                .count()
                .max(cluster_num - 1)
                + 1;
            debug!("Number of clusters:{}", cluster_num);
            let ((data, short), label, forbidden, answer) =
                select_within(region, data, label, forbidden, answer);
            debug!("Number of Reads:{}", data.len());
            assert_eq!(data.len(), forbidden.len());
            assert_eq!(data.len(), label.len() + answer.len());
            let predictions = {
                let (da, la, fo, an) = (&data, &label, &forbidden, &answer);
                // clustering(da, (la, an), fo, cluster_num, limit, c)
                let mut pred = clustering(da, (la, an), fo, cluster_num, limit, c);
                let pred_short = predict(&data, &pred, cluster_num, c, &short, 0);
                pred.extend(pred_short);
                pred
            };
            // assert_eq!(predictions.len(), data.len() + short.len());
            (0..cluster_num)
                .map(|cluster_idx| {
                    let cluster_idx = cluster_idx as u8;
                    predictions
                        .iter()
                        .zip(data.iter().chain(short.iter()))
                        // .zip(data.iter())
                        .filter(|(&p, _)| p == cluster_idx)
                        .map(|(_, read)| read.id().to_string())
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
    let mut components: Vec<HashSet<String>> = (0..(cluster_num * windows.len()))
        .filter_map(|parent| {
            if fu.find(parent).unwrap() != parent {
                return None;
            }
            debug!("find_cluseter");
            let component = clusterings
                .iter()
                .enumerate()
                .map(|(w, clusters)| {
                    clusters
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| {
                            let p = fu.find(j + w * cluster_num);
                            parent == p.unwrap_or(std::usize::MAX)
                        })
                        .fold(HashSet::new(), |mut acc, (j, cluster)| {
                            if cluster.len() > find_breakpoint::COVERAGE_THR {
                                info!("{}:{}", w, j);
                            }
                            acc.extend(cluster.clone());
                            acc
                        })
                })
                .fold(HashSet::new(), |mut acc, cluster| {
                    acc.extend(cluster);
                    acc
                });
            Some(component)
        })
        .filter(|component| component.len() > find_breakpoint::COVERAGE_THR)
        .collect();
    // And futher merging.
    'merge: loop {
        let len = components.len();
        debug!("Current Cluster:{}", len);
        for i in 0..len {
            for j in (i + 1)..len {
                for (k, initial_cluster) in initial_clusters.iter().enumerate() {
                    let reads = initial_cluster.ids();
                    let inter_i = components[i].intersection(&reads).count();
                    let inter_j = components[j].intersection(&reads).count();
                    debug!("{} shares {} reads with the init-cluster {}", i, inter_i, k);
                    debug!("{} shares {} reads with the init-cluster {}", j, inter_j, k);
                    if inter_i as f64 > CONNECTION_THR && inter_j as f64 > CONNECTION_THR {
                        let from = components.remove(j);
                        components[i].extend(from);
                        continue 'merge;
                    }
                }
            }
        }
        break;
    }
    debug!("Resulting in {} clusters.", components.len());
    let result: HashMap<String, u8> =
        components
            .into_iter()
            .enumerate()
            .fold(HashMap::new(), |mut acc, (idx, cluster)| {
                for read_name in cluster {
                    *acc.entry(read_name).or_default() = idx as u8;
                }
                acc
            });
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
    labels: (&[u8], &[u8]),
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    limit: u64,
    c: &Config,
) -> Vec<u8> {
    assert_eq!(forbidden.len(), data.len());
    assert_eq!((labels.0).len() + (labels.1).len(), data.len());
    use poa_clustering::{gibbs_sampling, DEFAULT_ALN};
    gibbs_sampling(data, labels, forbidden, cluster_num, limit, c, &DEFAULT_ALN)
}

pub fn predict(
    data: &[ERead],
    labels: &[u8],
    cluster_num: usize,
    c: &Config,
    input: &[ERead],
    seed: u64,
) -> Vec<u8> {
    assert_eq!(data.len(), labels.len());
    use poa_clustering::DEFAULT_ALN;
    poa_clustering::predict(data, labels, cluster_num, c, input, seed, &DEFAULT_ALN)
}
