#![allow(dead_code)]
#[macro_use]
extern crate log;
extern crate nalgebra as na;
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
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::{HashMap, HashSet};
pub mod annotate_contigs_to_reference;
pub mod assemble;
pub mod d3_data;
pub mod error_profile;
pub mod poa_clustering;
use poa_hmm::Config;
pub mod variant_calling;
// mod digamma;
// const WINDOW_SIZE: usize = 300;
// const OVERLAP: usize = 50;
const WINDOW_SIZE: usize = 20;
const OVERLAP: usize = 0;
const MIN_LEN: usize = 6_000;
const CONNECTION_THR: f64 = 0.5;
const MERGE_THR: usize = 50;
const NG_THR: usize = 10;
type Read<'a> = Vec<(usize, &'a [u8])>;
/// Main method. Decomposing the reads.
/// You should call "merge" method separatly(?) -- should be integrated with this function.
pub fn decompose(
    encoded_reads: Vec<last_tiling::EncodedRead>,
    initial_clusters: &[Cluster],
    contigs: &last_tiling::Contigs,
    config: &Config,
    cluster_num: usize,
    limit: u64,
) -> assemble::DecomposedResult {
    let ereads: Vec<_> = encoded_reads.iter().map(ERead::new_no_gapfill).collect();
    let coverages = get_coverages(contigs, &ereads);
    let (dataset, labels, forbidden) = initial_clustering(ereads, initial_clusters);
    let total_units = dataset.iter().map(|read| read.seq().len()).sum::<usize>();
    debug!("{} reads and {} units.", dataset.len(), total_units);
    debug!(
        "Forbids:{}",
        forbidden.values().filter(|e| !e.is_empty()).count()
    );
    let contigs: Vec<_> = (0..contigs.get_num_of_contigs())
        .map(|e| contigs.get_last_unit(e as u16).unwrap() as usize + 1)
        .collect();
    let windows: Vec<_> = contigs
        .iter()
        .zip(coverages)
        .enumerate()
        .flat_map(|(idx, (len, cov))| create_windows(idx, *len, &cov))
        .collect();
    let predicts = clustering_chunking(
        &dataset,
        &labels,
        &forbidden,
        cluster_num,
        config,
        limit,
        &windows,
    );
    let labels: HashMap<_, _> = dataset
        .iter()
        .zip(labels)
        .map(|(r, l)| (r.id.to_string(), l))
        .collect();
    let mut chunked_reads: Vec<_> = encoded_reads
        .iter()
        .flat_map(|r| {
            let label = labels.get(&r.id);
            let forbs = forbidden.get(&r.id);
            // There might be gappy reads.
            let entries = predicts.get(&r.id)?;
            Some(assemble::ChunkedRead::from(r, label, forbs, entries))
        })
        .collect();
    assemble::correct_reads::correct_reads(&mut chunked_reads);
    assemble::assemble_reads(chunked_reads)
}

fn initial_clustering(
    encoded_reads: Vec<ERead>,
    initial_clusters: &[Cluster],
) -> (Vec<ERead>, Vec<u8>, HashMap<String, Vec<u8>>) {
    let datasize = encoded_reads.len();
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
                .next();
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
    (dataset, labels, forbidden)
}

fn get_coverages(contigs: &last_tiling::Contigs, reads: &[ERead]) -> Vec<Vec<u32>> {
    let mut coverage: Vec<Vec<_>> = contigs
        .get_last_units()
        .into_iter()
        .map(|len| vec![0; len as usize + 1])
        .collect();
    for read in reads {
        for unit in read.seq() {
            let c = unit.contig() as usize;
            let u = unit.unit() as usize;
            coverage[c][u] += 1;
        }
    }
    coverage
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

pub fn create_windows(idx: usize, len: usize, covs: &[u32]) -> Vec<(u16, u16, u16)> {
    let mean = covs.iter().sum::<u32>() / covs.len() as u32;
    let mean_sq = covs.iter().fold(0, |x, y| x + y * y) / covs.len() as u32;
    let sd = ((mean_sq - mean * mean) as f64).sqrt().floor() as u32;
    //let thr = mean.max(3 * sd) - 3 * sd;
    let thr = mean / 4;
    debug!("Mean,SD,Thr={},{},{}", mean, sd, thr);
    let sub_windows: Vec<_> = {
        let mut sub_windows = vec![];
        let mut start = 0;
        while start < len {
            let end = start + covs[start..].iter().take_while(|&&c| c > thr).count();
            sub_windows.push((start, end));
            start = end + 1;
        }

        sub_windows
            .into_iter()
            .filter(|(s, e)| e - s > WINDOW_SIZE / 3)
            .collect()
    };
    sub_windows
        .into_iter()
        .flat_map(|(start, end)| {
            if end - start < WINDOW_SIZE {
                return vec![(start, end)];
            }
            let window_num = (end - start) / (WINDOW_SIZE - OVERLAP);
            let last_pos = start + (window_num - 1) * (WINDOW_SIZE - OVERLAP);
            assert!(last_pos >= 1);
            if end - last_pos < WINDOW_SIZE / 2 {
                (0..window_num)
                    .map(|i| {
                        let s = start + i * (WINDOW_SIZE - OVERLAP);
                        if i < last_pos - 1 {
                            (s, s + WINDOW_SIZE)
                        } else {
                            (s, end)
                        }
                    })
                    .collect::<Vec<_>>()
            } else {
                (0..=window_num)
                    .map(|i| {
                        let s = start + i * (WINDOW_SIZE - OVERLAP);
                        if i < last_pos {
                            (s, s + WINDOW_SIZE)
                        } else {
                            (s, end)
                        }
                    })
                    .collect::<Vec<_>>()
            }
        })
        .map(|(s, e)| (idx as u16, s as u16, e as u16))
        .collect()
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
            let mut forbid = forbidden.get(read.id()).cloned().unwrap_or_else(Vec::new);
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

#[derive(Clone, Debug)]
pub struct Entry<'a> {
    pub id: &'a str,
    pub window: usize,
    pub window_range: (u16, u16, u16),
    pub read_position: usize,
    pub seq: Vec<(usize, &'a [u8])>,
    pub forbid: &'a [u8],
    pub label: Option<u8>,
    pub assignment: u8,
}
impl<'a> Entry<'a> {
    fn new(
        id: &'a str,
        window: usize,
        window_range: (u16, u16, u16),
        read_position: usize,
        seq: Vec<(usize, &'a [u8])>,
        forbid: &'a [u8],
        label: Option<u8>,
    ) -> Self {
        Self {
            id,
            window,
            window_range,
            read_position,
            seq,
            forbid,
            label,
            assignment: 0,
        }
    }
}

pub fn clustering_chunking<'a>(
    data: &'a [ERead],
    label: &[u8],
    forbidden: &'a HashMap<String, Vec<u8>>,
    cluster_num: usize,
    c: &Config,
    limit: u64,
    windows: &[(u16, u16, u16)],
) -> HashMap<String, Vec<Entry<'a>>> {
    let mut pileups: Vec<Vec<_>> = vec![vec![]; windows.len()];
    for (pos, &(contig, start, end)) in windows.iter().enumerate() {
        for (idx, read) in data.iter().take(4000).enumerate() {
            let contained_window = read
                .seq
                .iter()
                .any(|u| u.contig == contig && start <= u.unit && u.unit < end);
            if contained_window {
                let filtered = read
                    .seq
                    .iter()
                    .enumerate()
                    .filter(|(_, u)| u.contig == contig && start <= u.unit && u.unit < end);
                let read_pos: usize = filtered.clone().next().unwrap().0;
                let seq: Vec<_> = filtered
                    .clone()
                    .map(|(_, u)| ((u.unit - start) as usize, u.bases.as_slice()))
                    .collect();
                let lab = label.get(idx).cloned();
                let forb = forbidden.get(&read.id).unwrap();
                let range = (contig, start, end);
                let entry = Entry::new(&read.id, pos, range, read_pos, seq, forb, lab);
                pileups[pos].push(entry);
            }
        }
    }
    // Put segments with labels forward.
    pileups.par_iter_mut().for_each(|pileup| {
        use std::cmp::Ordering::*;
        pileup.sort_by(|a, b| match (a.label.is_some(), b.label.is_some()) {
            (true, true) | (false, false) => Equal,
            (true, false) => Less,
            (false, true) => Greater,
        });
    });
    // Parallelize here to get most efficient algorithm.
    pileups
        .par_iter_mut()
        .zip(windows.into_par_iter())
        .for_each(|(pileup, range)| {
            let label_map = pileup.iter().fold(HashMap::new(), |mut res, entry| {
                let next = res.len() as u8;
                if let Some(cluster) = entry.label {
                    res.entry(cluster).or_insert(next);
                }
                res
            });
            // This does not break correspondance between entries and labels.
            let labels: Vec<_> = pileup
                .iter()
                .filter_map(|entry| entry.label)
                .map(|cl| label_map[&cl])
                .collect();
            debug!("Start {:?}", range);
            let coverage = pileup.len();
            debug!("{} labels. Coverage is {}", label.len(), coverage);
            let forbs: Vec<Vec<u8>> = pileup
                .iter()
                .map(|entry| {
                    entry
                        .forbid
                        .iter()
                        .filter_map(|x| label_map.get(x))
                        .copied()
                        .collect()
                })
                .collect();
            let cluster_num = label_map.len().max(cluster_num);
            let chain_len = (range.2 - range.1) as usize;
            let data: Vec<_> = pileup.iter().map(|e| e.seq.clone()).collect();
            let predictions = poa_clustering::gibbs_sampling(
                &data,
                &labels,
                None,
                &forbs,
                chain_len,
                cluster_num,
                limit,
                c,
                &poa_clustering::DEFAULT_ALN,
                coverage,
            );
            pileup
                .iter_mut()
                .zip(predictions)
                .for_each(|(mut e, p)| e.assignment = p);
        });
    if log_enabled!(log::Level::Debug) {
        let id2desc: HashMap<_, _> = data
            .iter()
            .filter_map(|r| r.desc.as_ref().map(|desc| (r.id.to_string(), desc)))
            .collect();
        for (idx, pileup) in pileups.iter().enumerate() {
            for entry in pileup.iter() {
                let desc = match id2desc.get(entry.id) {
                    Some(res) => res,
                    None => continue,
                };
                debug!("{}\t{}\t{}", idx, entry.assignment, desc);
            }
        }
    }
    let mut encoded_reads: HashMap<_, Vec<_>> = HashMap::new();
    for pileup in pileups {
        for entry in pileup {
            encoded_reads
                .entry(entry.id.to_string())
                .or_default()
                .push(entry);
        }
    }
    encoded_reads
        .iter_mut()
        .for_each(|(_, entries)| entries.sort_by_key(|e| e.read_position));
    encoded_reads
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
    _forbidden: &HashMap<String, Vec<u8>>,
) -> Vec<HashSet<String>> {
    let (components, background): (Vec<_>, Vec<_>) =
        components.into_iter().partition(|component| {
            initial_clusters
                .iter()
                .any(|cl| cl.ids().intersection(component).count() > MERGE_THR)
        });
    let mut background: HashSet<String> =
        background.into_iter().fold(HashSet::new(), |mut x, y| {
            x.extend(y);
            x
        });
    let mut result = vec![];
    let is_occupying = determine_occupy(&components, data);
    for (component, is_occupying) in components.into_iter().zip(is_occupying) {
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
    let component_len = components.len();
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
    let mean = total.iter().sum::<usize>() / total.len();
    let avesq = total.iter().map(|e| e * e).sum::<usize>() / total.len();
    let thr = mean.max(avesq * 3) - (avesq * 3);
    (0..component_len)
        .map(|cl| {
            count
                .iter()
                .zip(total.iter())
                .filter(|x| *x.1 > thr)
                .any(|(pos, &tot)| pos[cl] > 8 * tot / 10)
        })
        .collect()
}

// fn resume_clustering(
//     clusterings: Vec<Vec<HashSet<String>>>,
//     windowlen: usize,
//     forbidden: &HashMap<String, Vec<u8>>,
//     initial_clusters: &[Cluster],
//     data: &[ERead],
// ) -> Vec<Option<u8>> {
//     debug!(
//         "On windows:{}",
//         clusterings.iter().map(|cl| cl.len()).sum::<usize>()
//     );
//     let components =
//         merge_windows_by_bipartite(clusterings, windowlen, forbidden, initial_clusters, data);
//     debug!("Merged:{}", components.len());
//     let components = merge_by_initial_cluster(forbidden, initial_clusters, components);
//     debug!("Merge by Init:{}", components.len());
//     let components = if initial_clusters.is_empty() {
//         components
//     } else {
//         merge_with_background(data, initial_clusters, components, forbidden)
//     };
//     debug!("Resulting in {} clusters.", components.len());
//     let result: HashMap<String, u8> =
//         components
//             .into_iter()
//             .enumerate()
//             .fold(HashMap::new(), |mut acc, (idx, cluster)| {
//                 for read_name in cluster {
//                     acc.insert(read_name, idx as u8);
//                 }
//                 acc
//             });
//     data.iter().map(|r| result.get(r.id()).cloned()).collect()
// }

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
    // debug!("Shares {} reads, NGS {} reads.", share, ngs);
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

// Return clustering. If the value is None,
// the element can not be assigned to any cluster.
pub fn clustering(
    _data: &[ERead],
    _labels: &[u8],
    _forbidden: &[Vec<u8>],
    _cluster_num: usize,
    _limit: u64,
    _c: &Config,
    _coverage: usize,
) -> Vec<Option<u8>> {
    vec![]
    // assert_eq!(forbidden.len(), data.len());
    // use poa_clustering::{gibbs_sampling, DEFAULT_ALN};
    // let answer = vec![0u8; data.len() - labels.len()];
    // let (ls, cl) = ((labels, answer.as_slice()), cluster_num);
    // gibbs_sampling(data, ls, forbidden, cl, limit, c, &DEFAULT_ALN, coverage)
}
