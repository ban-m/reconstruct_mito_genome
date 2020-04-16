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
use poa_hmm::Config;
pub mod variant_calling;
const WINDOW_SIZE: usize = 200;
const OVERLAP: usize = 25;
const MIN_LEN: usize = 5_000;
const CONNECTION_THR: f64 = 0.7;
const MERGE_THR: usize = 5;
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
    cluster_num: usize,
    limit: u64,
) -> Vec<(String, Option<u8>)> {
    let (mut dataset, lab_ans, masked_region, forbidden) =
        initial_clustering(encoded_reads, initial_clusters, contigs, repeats);
    let (labels, answer) = lab_ans;
    let total_units = dataset.iter().map(|read| read.seq().len()).sum::<usize>();
    debug!("{} reads and {} units.", dataset.len(), total_units);
    let forbids = forbidden.values().filter(|e| !e.is_empty()).count();
    debug!("Forbids:{}", forbids);
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
    // Remove chunks from masked region.
    dataset.iter_mut().for_each(|read| {
        read.seq = read
            .seq
            .iter()
            .filter(|u| !masked_region[u.contig()][u.unit()])
            .cloned()
            .collect();
    });
    let predicts = clustering_chunking(
        &dataset,
        (&labels, &answer),
        &forbidden,
        initial_clusters,
        cluster_num,
        &contigs,
        config,
        limit,
        &windows,
    );
    dataset
        .into_iter()
        .zip(predicts.iter())
        .map(|(read, cl)| (read.id().to_string(), *cl))
        .collect()
}

pub fn resume_decompose(
    encoded_reads: Vec<ERead>,
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
    resume: Vec<Vec<HashSet<String>>>,
) -> Vec<(String, Option<u8>)> {
    let (dataset, _, masked_region, forbidden) =
        initial_clustering(encoded_reads, initial_clusters, contigs, repeats);
    let windowlen = (0..contigs.get_num_of_contigs())
        .map(|e| contigs.get_last_unit(e as u16).unwrap() as usize + 1)
        .zip(masked_region.iter())
        .enumerate()
        .flat_map(|(idx, (len, mask))| create_windows(idx, len, mask))
        .count();
    let predicts = resume_clustering(resume, windowlen, &forbidden, initial_clusters, &dataset);
    dataset
        .into_iter()
        .zip(predicts.iter())
        .map(|(read, cl)| (read.id().to_string(), *cl))
        .collect()
}

fn initial_clustering(
    encoded_reads: Vec<ERead>,
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
) -> (
    Vec<ERead>,
    (Vec<u8>, Vec<u8>),
    Vec<Vec<bool>>,
    HashMap<String, Vec<u8>>,
) {
    let datasize = encoded_reads.len();
    let masked_region = get_masked_region(&initial_clusters, &contigs, repeats);
    let mut forbidden: HashMap<String, _> = HashMap::new();
    let mut labels: Vec<_> = vec![];
    let (assigned_reads, unassigned_reads): (Vec<_>, Vec<_>) = encoded_reads
        .into_iter()
        .filter_map(|mut read| {
            let forbid: Vec<_> = initial_clusters
                .iter()
                .filter(|cluster| cluster.is_spanned_by(&read))
                .map(|c| c.id as u8)
                .collect();
            let seq: Vec<_> = read.seq().iter().cloned().collect();
            forbidden.insert(read.id().to_string(), forbid);
            if seq.is_empty() {
                None
            } else {
                *read.seq_mut() = seq;
                Some(read)
            }
        })
        .partition(|read| {
            let matched_cluster = initial_clusters
                .iter()
                .filter_map(|cr| if cr.has(read.id()) { Some(cr.id) } else { None })
                .nth(0);
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
    let answer = vec![0; unassigned_reads.len()];
    let dataset = vec![assigned_reads, unassigned_reads].concat();
    (dataset, (labels, answer), masked_region, forbidden)
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
            start = end - OVERLAP;
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
    forbidden: &HashMap<String, Vec<u8>>,
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
            let forbid = forbidden.get(read.id()).cloned().unwrap_or(vec![]);
            s_forbid.push(forbid);
            s_data.push(read);
            if idx < border {
                s_label.push(label[idx]);
            } else {
                s_answer.push(answer[idx - border]);
            }
        } else if read.seq().len() > original_len / 2 {
            shorter_reads.push(read);
        }
    }
    assert_eq!(s_data.len(), s_forbid.len());
    ((s_data, shorter_reads), s_label, s_forbid, s_answer)
}

fn find_matching(
    prev: &[HashSet<String>],
    after: &[HashSet<String>],
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    lengths: &HashMap<String, usize>,
) -> Vec<(usize, usize)> {
    let node1 = prev.len();
    let node2 = after.len();
    let filter_short = |hm: &HashSet<String>| {
        hm.iter()
            .filter(|&id| match lengths.get(id) {
                Some(&res) => res > MIN_LEN / last_tiling::UNIT_SIZE,
                None => false,
            })
            .cloned()
            .collect::<HashSet<_>>()
    };
    let prev: Vec<HashSet<String>> = prev.iter().map(filter_short).collect();
    let after: Vec<HashSet<String>> = after.iter().map(filter_short).collect();
    debug!("Dump graph");
    let prev_boundary: Vec<HashSet<_>> = prev
        .iter()
        .map(|from| {
            after.iter().fold(HashSet::new(), |mut acc, x| {
                acc.extend(x.intersection(from).cloned());
                acc
            })
        })
        .collect();
    let after_boundary: Vec<HashSet<_>> = after
        .iter()
        .map(|to| {
            prev.iter().fold(HashSet::new(), |mut acc, x| {
                acc.extend(x.intersection(to).cloned());
                acc
            })
        })
        .collect();
    let graph: Vec<Vec<(usize, f64)>> = prev
        .iter()
        .zip(prev_boundary)
        .enumerate()
        .map(|(from, (cl1, cl1_boundary))| {
            after
                .iter()
                .zip(after_boundary.iter())
                .enumerate()
                .filter(|&(to, (cl2, _))| {
                    is_mergiable(from, to, cl1, cl2, forbidden, initial_clusters)
                })
                .filter_map(|(to, (cl2, cl2_boundary))| {
                    let intersect = cl1.intersection(cl2).count();
                    let union = cl1_boundary.union(cl2_boundary).count();
                    let sim = intersect as f64 / union as f64;
                    debug!("{}->({:.3}={}/{})->{}", from, sim, intersect, union, to);
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
    debug!("Path Selected.");
    for &(from, to) in &res {
        debug!("{}-{}", from, to);
    }
    res
}

fn is_mergiable(
    from: usize,
    to: usize,
    cl1: &HashSet<String>,
    cl2: &HashSet<String>,
    forbidden: &HashMap<String, Vec<u8>>,
    clusters: &[Cluster],
) -> bool {
    let cl1_cluster: HashSet<_> = clusters
        .iter()
        .filter_map(|cl| {
            if cl1.iter().any(|id| cl.has(id)) {
                Some(cl.id as u8)
            } else {
                None
            }
        })
        .collect();
    let cl2_cluster: HashSet<_> = clusters
        .iter()
        .filter_map(|cl| {
            if cl2.iter().any(|id| cl.has(id)) {
                Some(cl.id as u8)
            } else {
                None
            }
        })
        .collect();
    // Check whether cl1_cluster violates cl2's forbidden set
    // or cl2_cluster violates cl1's forbidden set.
    let cl1_violate: usize = cl1
        .iter()
        .filter_map(|id| forbidden.get(id))
        .filter(|forbid| forbid.iter().any(|f| cl2_cluster.contains(f)))
        .count();
    let cl2_violate: usize = cl2
        .iter()
        .filter_map(|id| forbidden.get(id))
        .filter(|forbid| forbid.iter().any(|f| cl1_cluster.contains(f)))
        .count();
    if cl1_violate < 10 && cl2_violate < 10 {
        true
    } else {
        debug!("{}({:?})-{}({:?})", from, cl1_cluster, to, cl2_cluster);
        debug!("Violates:{}/{}", cl1_violate, cl2_violate);
        false
    }
}

fn map_labels(label: &[u8]) -> HashMap<u8, u8> {
    let mut clusters: Vec<_> = label.iter().copied().collect();
    clusters.sort();
    clusters.dedup();
    clusters
        .into_iter()
        .enumerate()
        .map(|(to, from)| (from, to as u8))
        .collect()
}

/// Clustering after chunking the reference into several chunks.
pub fn clustering_chunking(
    data: &[ERead],
    (label, answer): (&[u8], &[u8]),
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
    limit: u64,
    windows: &[(u16, u16, u16)],
) -> Vec<Option<u8>> {
    debug!("The contig lengths:{:?}", contigs);
    debug!("There are {} windows in total.", windows.len());
    for w in windows {
        debug!("{:?}", w);
    }
    let clusterings: Vec<Vec<HashSet<String>>> = windows
        .iter()
        .enumerate()
        .map(|(idx, &region)| {
            let (contig, start, end) = region;
            debug!("{}/{}:{}-{}", idx, contig, start, end);
            let ((data, short), label, forbidden, answer) =
                select_within(region, data, label, forbidden, answer);
            debug!("Number of Reads:{}", data.len());
            assert_eq!(data.len(), forbidden.len());
            assert_eq!(data.len(), label.len() + answer.len());
            let cluster_map: HashMap<u8, u8> = map_labels(&label);
            let cluster_num = cluster_map.len().max(cluster_num - 1) + 1;
            debug!("Number of clusters:{}", cluster_num);
            debug!("Mapping:{:?}", cluster_map);
            let label: Vec<_> = label.iter().map(|c| cluster_map[c]).collect();
            let forbidden: Vec<Vec<_>> = forbidden
                .iter()
                .map(|f| {
                    f.iter()
                        .filter_map(|c| cluster_map.get(c))
                        .copied()
                        .collect()
                })
                .collect();
            let forbs = forbidden.iter().filter(|e| !e.is_empty()).count();
            debug!("{} labels and {} noempty forbidden", label.len(), forbs);
            let predictions = {
                let (da, la, fo, an) = (&data, &label, &forbidden, &answer);
                let mut pred = clustering(da, (la, an), fo, cluster_num, limit, c);
                let pred_short = predict(da, &pred, cluster_num, c, &short, 0);
                pred.extend(pred_short);
                pred
            };
            dump_pred(&predictions, &data, &short, idx);
            (0..cluster_num)
                .map(|cluster_idx| {
                    let cluster_idx = cluster_idx as u8;
                    predictions
                        .iter()
                        .zip(data.iter().chain(short.iter()))
                        .filter_map(|(p, x)| p.map(|p| (p, x)))
                        .filter(|&(p, _)| p == cluster_idx)
                        .map(|(_, read)| read.id().to_string())
                        .collect()
                })
                .collect()
        })
        .collect();
    resume_clustering(
        clusterings,
        windows.len(),
        forbidden,
        initial_clusters,
        data,
    )
}

fn resume_clustering(
    clusterings: Vec<Vec<HashSet<String>>>,
    windowlen: usize,
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    data: &[ERead],
) -> Vec<Option<u8>> {
    let max_cluster_num = clusterings.iter().map(|e| e.len()).max().unwrap_or(0);
    let lengths: HashMap<String, usize> = data
        .iter()
        .map(|r| (r.id().to_string(), r.seq.len()))
        .collect();
    let mut fu = find_union::FindUnion::new(max_cluster_num * windowlen);
    for idx in 0..clusterings.len() {
        let prev_idx = idx;
        let after_idx = (idx + 1) % clusterings.len();
        let prev = &clusterings[prev_idx];
        let after = &clusterings[after_idx];
        for (i, j) in find_matching(prev, after, forbidden, initial_clusters, &lengths) {
            let i = i + prev_idx * max_cluster_num;
            let j = j + after_idx * max_cluster_num;
            fu.unite(i, j).unwrap();
        }
    }
    // Then, iteratively take components.
    let mut components: Vec<HashSet<String>> = (0..(max_cluster_num * windowlen))
        .filter_map(|parent| {
            if fu.find(parent).unwrap() != parent {
                return None;
            }
            let component = clusterings
                .iter()
                .enumerate()
                .map(|(w, clusters)| {
                    clusters
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| {
                            let p = fu.find(j + w * max_cluster_num);
                            parent == p.unwrap_or(std::usize::MAX)
                        })
                        .fold(HashSet::new(), |mut acc, (_, cluster)| {
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
                    if inter_i > MERGE_THR && inter_j > MERGE_THR {
                        debug!("{} shares {} reads with the init-cluster {}", i, inter_i, k);
                        debug!("{} shares {} reads with the init-cluster {}", j, inter_j, k);
                        let from = components.remove(j);
                        components[i].extend(from);
                        continue 'merge;
                    }
                }
            }
        }
        break;
    }
    // Filtering out clusters which do not overlap any initial clusters.
    let components: Vec<HashSet<String>> = if initial_clusters.is_empty() {
        components
    } else {
        components
            .into_iter()
            .filter(|component| {
                initial_clusters
                    .iter()
                    .any(|cl| cl.ids().intersection(component).count() > 10)
            })
            .collect()
    };
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
    data.iter().map(|r| result.get(r.id()).cloned()).collect()
}

fn dump_pred(assignments: &[Option<u8>], data: &[ERead], data2: &[ERead], idx: usize) {
    for (asn, data) in assignments.iter().zip(data.iter().chain(data2.iter())) {
        let no_desc = String::from("None");
        let desc = data.desc().unwrap_or(&no_desc);
        let id = data.id();
        match asn {
            Some(res) => debug!("PRED\t{}\t{}\t{}\t{}", idx, res, id, desc),
            None => debug!("PRED\t{}\tNone\t{}\t{}", idx, id, desc),
        }
    }
}

/// Return clustering. If the value is None,
/// the element can not be assigned to any cluster.
pub fn clustering(
    data: &[ERead],
    labels: (&[u8], &[u8]),
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    limit: u64,
    c: &Config,
) -> Vec<Option<u8>> {
    assert_eq!(forbidden.len(), data.len());
    assert_eq!((labels.0).len() + (labels.1).len(), data.len());
    use poa_clustering::{gibbs_sampling, DEFAULT_ALN};
    gibbs_sampling(data, labels, forbidden, cluster_num, limit, c, &DEFAULT_ALN)
}

pub fn predict(
    data: &[ERead],
    labels: &[Option<u8>],
    cluster_num: usize,
    c: &Config,
    input: &[ERead],
    seed: u64,
) -> Vec<Option<u8>> {
    assert_eq!(data.len(), labels.len());
    use poa_clustering::DEFAULT_ALN;
    poa_clustering::predict(data, labels, cluster_num, c, input, seed, &DEFAULT_ALN)
}
