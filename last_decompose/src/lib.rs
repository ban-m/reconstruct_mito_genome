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
// mod digamma;
const WINDOW_SIZE: usize = 250;
const OVERLAP: usize = 60;
const MIN_LEN: usize = 6_000;
const CONNECTION_THR: f64 = 0.8;
const MERGE_THR: usize = 50;
const NG_THR: usize = 10;
type Read = Vec<(usize, Vec<u8>)>;
/// Main method. Decomposing the reads.
/// You should call "merge" method separatly(?) -- should be integrated with this function.
pub fn decompose(
    encoded_reads: Vec<ERead>,
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
    config: &Config,
    cluster_num: usize,
    limit: u64,
) -> Vec<(String, Option<u8>)> {
    let (dataset, labels, masked_region, forbidden) =
        initial_clustering(encoded_reads, initial_clusters, contigs, repeats);
    let total_units = dataset.iter().map(|read| read.seq().len()).sum::<usize>();
    debug!("{} reads and {} units.", dataset.len(), total_units);
    let forbids = forbidden.values().filter(|e| !e.is_empty()).count();
    debug!("Forbids:{}", forbids);
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
    // dataset.iter_mut().for_each(|read| {
    //     read.seq = read
    //         .seq
    //         .iter()
    //         .filter(|u| !masked_region[u.contig()][u.unit()])
    //         .cloned()
    //         .collect();
    // });
    let predicts = clustering_chunking(
        &dataset,
        &labels,
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
    Vec<u8>,
    Vec<Vec<bool>>,
    HashMap<String, Vec<u8>>,
) {
    let datasize = encoded_reads.len();
    let masked_region = get_masked_region(&initial_clusters, &contigs, repeats);
    let mut forbidden: HashMap<String, _> = HashMap::new();
    let mut labels: Vec<_> = vec![];
    let (assigned_reads, unassigned_reads): (Vec<_>, Vec<_>) = encoded_reads
        .into_iter()
        .inspect(|read| {
            let forbid: Vec<_> = initial_clusters
                .iter()
                .filter(|cluster| cluster.is_spanned_by(&read))
                .map(|c| c.id as u8)
                .collect();
            forbidden.insert(read.id().to_string(), forbid);
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
    assert_eq!(assigned_reads.len() + unassigned_reads.len(), datasize);
    debug!(
        "Unassigned reads before removing contained:{}",
        unassigned_reads.len()
    );
    let dataset = vec![assigned_reads, unassigned_reads].concat();
    (dataset, labels, masked_region, forbidden)
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
    init_cluster: &[Cluster],
) -> (Vec<ERead>, Vec<u8>, Vec<Vec<u8>>) {
    let (mut s_data, mut s_label, mut s_forbid) = (vec![], vec![], vec![]);
    debug!("Selecting {}\t{}\t{}...", contig, start, end);
    let additional_forbiddens = {
        let range = (contig, start, end);
        let clusters: Vec<_> = init_cluster.iter().filter(|cr| cr.overlap(range)).collect();
        move |read: &ERead| -> Vec<u8> {
            clusters
                .iter()
                .filter(|cr| cr.locally_forbids(range, read))
                .map(|cr| cr.id as u8)
                .collect()
        }
    };
    let border = label.len();
    for (idx, read) in data.iter().enumerate() {
        let original_len = read.seq().len();
        let read = read.clone_within(contig, start, end);
        let unit_thr = (MIN_LEN / last_tiling::UNIT_SIZE).min(original_len / 2);
        if read.seq().len() > unit_thr {
            let mut forbid = forbidden.get(read.id()).cloned().unwrap_or(vec![]);
            forbid.extend(additional_forbiddens(&read));
            s_forbid.push(forbid);
            s_data.push(read);
            if idx < border {
                s_label.push(label[idx]);
            }
        }
    }
    assert_eq!(s_data.len(), s_forbid.len());
    (s_data, s_label, s_forbid)
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
                    let intersect = cl1.intersection(cl2).count().pow(2);
                    let union = cl1_boundary.len() * cl2_boundary.len();
                    let sim = (intersect as f64 / union as f64).sqrt();
                    debug!("{}->({:.3}={}/{})->{}", from, sim, intersect, union, to);
                    if sim > CONNECTION_THR || union == 0 {
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
            if cl1.iter().filter(|id| cl.has(id)).count() > 10 {
                Some(cl.id as u8)
            } else {
                None
            }
        })
        .collect();
    let cl2_cluster: HashSet<_> = clusters
        .iter()
        .filter_map(|cl| {
            if cl2.iter().filter(|id| cl.has(id)).count() > 10 {
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
    label: &[u8],
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
    for (idx, w) in windows.iter().enumerate() {
        debug!("{}:{:?}", idx, w);
    }
    let clusterings: Vec<Vec<HashSet<String>>> = windows
        .iter()
        .enumerate()
        .map(|(idx, &region)| {
            let (contig, start, end) = region;
            debug!("{}/{}:{}-{}", idx, contig, start, end);
            let (data, label, forbidden) =
                select_within(region, data, label, forbidden, initial_clusters);
            debug!("Number of Reads:{}", data.len());
            assert_eq!(data.len(), forbidden.len());
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
            let len = (end - start) as usize;
            let coverage = data.iter().map(|r| r.seq.len()).sum::<usize>() / len;
            debug!("{} labels and {} noempty forbidden", label.len(), forbs);
            debug!("Coverage is {}", coverage);
            let predictions =
                clustering(&data, &label, &forbidden, cluster_num, limit, c, coverage);
            dump_pred(&predictions, &data, idx);
            (0..cluster_num)
                .map(|cluster_idx| {
                    let cluster_idx = cluster_idx as u8;
                    predictions
                        .iter()
                        .zip(data.iter())
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

fn merge_windows_by_bipartite(
    clusterings: Vec<Vec<HashSet<String>>>,
    windowlen: usize,
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    data: &[ERead],
) -> Vec<HashSet<String>> {
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
    (0..(max_cluster_num * windowlen))
        .filter_map(|parent| {
            if fu.find(parent).unwrap() == parent {
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
            } else {
                None
            }
        })
        .collect()
}

fn merge_by_initial_cluster(
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    mut components: Vec<HashSet<String>>,
) -> Vec<HashSet<String>> {
    for cluster in initial_clusters {
        debug!("Merging by {:?}", cluster);
        components = {
            let (merged, mut result): (Vec<_>, Vec<_>) = components
                .into_iter()
                .partition(|c| is_overlap(c, cluster, forbidden));
            debug!("There are {} cluster merged.", merged.len());
            let merged = merged.into_iter().fold(HashSet::new(), |mut acc, x| {
                acc.extend(x);
                acc
            });
            result.push(merged);
            result
        };
    }
    components
}

fn merge_with_background(
    data: &[ERead],
    initial_clusters: &[Cluster],
    components: Vec<HashSet<String>>,
    forbidden: &HashMap<String, Vec<u8>>,
) -> Vec<HashSet<String>> {
    let (components, background): (Vec<_>, Vec<_>) =
        components.into_iter().partition(|component| {
            initial_clusters
                .iter()
                .any(|cl| cl.ids().intersection(component).count() > MERGE_THR)
        });
    let mut background = background.into_iter().fold(HashSet::new(), |mut x, y| {
        x.extend(y);
        x
    });
    let _forbidden: HashSet<_> = background
        .iter()
        .filter_map(|id| forbidden.get(id))
        .flat_map(|t| t)
        .copied()
        .collect();
    let is_occupying = determine_occupy(&components, data);
    let mut result = vec![];
    for (component, is_occupying) in components.into_iter().zip(is_occupying) {
        let _cluster: HashSet<_> = initial_clusters
            .iter()
            .filter(|cl| component.iter().filter(|id| cl.has(id)).count() > 10)
            .map(|cl| cl.id as u8)
            .collect();
        if is_occupying {
            background.extend(component);
        } else {
            result.push(component);
        }
    }
    result.push(background);
    result
}

fn determine_occupy(components: &Vec<HashSet<String>>, data: &[ERead]) -> Vec<bool> {
    let (matrix_pos, chain_len) = poa_clustering::to_pos(data);
    let components: HashMap<_, usize> =
        components
            .iter()
            .enumerate()
            .fold(HashMap::new(), |mut x, (idx, cl)| {
                for id in cl {
                    x.insert(id.clone(), idx);
                }
                x
            });
    let mut count: Vec<_> = vec![vec![0; components.len()]; chain_len];
    let mut total = vec![0; chain_len];
    for read in data.iter() {
        if let Some(&cl) = components.get(read.id()) {
            for unit in read.seq() {
                let pos = matrix_pos[unit.contig()][unit.unit()];
                count[pos][cl] += 1;
                total[pos] += 1;
            }
        } else {
            for unit in read.seq() {
                let pos = matrix_pos[unit.contig()][unit.unit()];
                total[pos] += 1;
            }
        }
    }
    let thr = total.iter().sum::<usize>() / total.len();
    (0..components.len())
        .map(|cl| {
            count
                .iter()
                .zip(total.iter())
                .filter(|x| *x.1 > thr)
                .any(|(pos, &tot)| pos[cl] > 8 * tot / 10)
        })
        .collect()
}

fn resume_clustering(
    clusterings: Vec<Vec<HashSet<String>>>,
    windowlen: usize,
    forbidden: &HashMap<String, Vec<u8>>,
    initial_clusters: &[Cluster],
    data: &[ERead],
) -> Vec<Option<u8>> {
    let components =
        merge_windows_by_bipartite(clusterings, windowlen, forbidden, initial_clusters, data);
    let components = merge_by_initial_cluster(forbidden, initial_clusters, components);
    let components = if initial_clusters.is_empty() {
        components
    } else {
        merge_with_background(data, initial_clusters, components, forbidden)
    };
    debug!("Resulting in {} clusters.", components.len());
    let result: HashMap<String, u8> =
        components
            .into_iter()
            .enumerate()
            .fold(HashMap::new(), |mut acc, (idx, cluster)| {
                for read_name in cluster {
                    acc.insert(read_name, idx as u8);
                }
                acc
            });
    data.iter().map(|r| result.get(r.id()).cloned()).collect()
}

fn is_overlap(
    component: &HashSet<String>,
    cluster: &Cluster,
    forbids: &HashMap<String, Vec<u8>>,
) -> bool {
    let ngs: usize = component
        .iter()
        .filter_map(|id| forbids.get(id))
        .filter(|forbid| forbid.contains(&(cluster.id as u8)))
        .count();
    let share = component.iter().filter(|id| cluster.has(id)).count();
    debug!("Shares {} reads, NGS {} reads.", share, ngs);
    share > MERGE_THR && ngs <= NG_THR
}

fn dump_pred(assignments: &[Option<u8>], data: &[ERead], idx: usize) {
    for (asn, data) in assignments.iter().zip(data.iter()) {
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
    labels: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    limit: u64,
    c: &Config,
    coverage: usize,
) -> Vec<Option<u8>> {
    assert_eq!(forbidden.len(), data.len());
    use poa_clustering::{gibbs_sampling, DEFAULT_ALN};
    let answer = vec![0u8; data.len() - labels.len()];
    let (ls, cl) = ((labels, answer.as_slice()), cluster_num);
    gibbs_sampling(data, ls, forbidden, cl, limit, c, &DEFAULT_ALN, coverage)
}
