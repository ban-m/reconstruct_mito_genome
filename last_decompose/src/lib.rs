#[macro_use]
extern crate log;
extern crate bio_utils;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate last_tiling;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
extern crate serde;
pub use find_breakpoint::critical_regions;
use rayon::prelude::*;
pub mod bipartite_matching;
mod eread;
pub mod find_breakpoint;
mod find_union;
pub mod utils;
use dbg_hmm::*;
pub use eread::*;
pub use find_breakpoint::CriticalRegion;
use find_breakpoint::ReadClassify;
use last_tiling::repeat::RepeatPairs;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand::{thread_rng, Rng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::{HashMap, HashSet};
pub mod error_profile;
const PRIOR: bool = true;
const MODEL_CHECK: bool = true;
const NUM_OF_BALL: usize = 100;
const WINDOW_SIZE: usize = 300;
const OVERLAP: usize = 100;
const LK_LIMIT: f64 = -140.;
// const LK_LIMIT: f64 = -60.;
const OFFSET: f64 = 0.;
const A_PRIOR: f64 = -0.23;
const B_PRIOR: f64 = 3.17;
const INIT_BETA: f64 = 0.008;
// const MAX_BETA: f64 = 1.;
// This is the search factor for the doubling-step.
const FACTOR: f64 = 1.4;
// Sampling times for Gibbs sampling.
const SAMPLING: usize = 90;
// This is the factor we multiply at each iteration for beta.
// Note that this factor also scaled for the maximum coverage.
const BETA_STEP: f64 = 1.1;
// Maximum beta. Until this beta, we multiply beta_step for the current beta.
// Loop number for Gibbs sampling.
const LOOP_NUM: usize = 20;
// Initial picking probability.
const INIT_PICK_PROB: f64 = 0.05;
const SOE_PER_DATA_ENTROPY: f64 = 0.20;
const K: usize = 6;
/// Main method. Decomposing the reads.
/// You should call "merge" method separatly(?) -- should be integrated with this function.
// TODO: Make a procedure to remove chimeric reads.
// -> If a read has 'junction', and in the "unassigned" bucket,
// it means such read can be safely removed.
pub fn decompose(
    encoded_reads: Vec<ERead>,
    critical_regions: &[CriticalRegion],
    contigs: &last_tiling::Contigs,
    repeats: &[RepeatPairs],
    config: dbg_hmm::Config,
    answer: &HashMap<String, u8>,
) -> Vec<(String, u8)> {
    for (idx, cr) in critical_regions.iter().enumerate() {
        debug!("{}\t{}", idx, cr);
    }
    let datasize = encoded_reads.len();
    let mut unassigned_reads: Vec<_> = vec![];
    let mut assigned_reads: Vec<_> = vec![];
    let mut labels: Vec<_> = vec![];
    for read in encoded_reads {
        let matched_cluster = critical_regions
            .iter()
            .enumerate()
            .filter(|(_, cr)| cr.along_with(&read))
            .nth(0);
        if let Some((idx, _)) = matched_cluster {
            assigned_reads.push(read);
            labels.push(idx as u8);
        } else {
            unassigned_reads.push(read);
        }
    }
    assert_eq!(labels.len(), assigned_reads.len());
    assert_eq!(assigned_reads.len() + unassigned_reads.len(), datasize);
    let masked_region = get_masked_region(&critical_regions, &contigs, repeats);
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
            let crs: Vec<_> = critical_regions
                .iter()
                .enumerate()
                .filter_map(|(idx, cr)| {
                    if cr.is_spanned_by(&read) {
                        Some(idx as u8)
                    } else {
                        None
                    }
                })
                .collect();
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
                forbidden.push(crs);
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
        K,
        critical_regions.len().max(1) + 1,
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
    critical_regions: &[CriticalRegion],
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
    let ranges: Vec<_> = critical_regions
        .iter()
        .flat_map(|cr| match cr {
            CriticalRegion::CP(ref cp) => vec![
                (cp.contig1().contig(), cp.contig1().range()),
                (cp.contig1().contig(), cp.contig2().range()),
            ],
            CriticalRegion::CR(ref cr) => vec![(cr.contig().contig(), cr.contig().range())],
        })
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

pub struct ModelFactory<'a> {
    // Contig -> Unit -> Seqs
    chunks: Vec<Vec<Vec<&'a [u8]>>>,
    weights: Vec<Vec<Vec<f64>>>,
    // Contig -> Unit
    factories: Vec<Vec<Factory>>,
    buffers: Vec<Vec<Vec<f64>>>,
    k: usize,
}

impl<'a> ModelFactory<'a> {
    pub fn new(contigs: &[usize], data: &'a [ERead], k: usize) -> Self {
        let mut chunks: Vec<Vec<Vec<&[u8]>>> = contigs.iter().map(|&e| vec![vec![]; e]).collect();
        let weights: Vec<Vec<Vec<f64>>> = contigs.iter().map(|&e| vec![vec![]; e]).collect();
        for read in data.iter() {
            for chunk in read.seq.iter() {
                chunks[chunk.contig()][chunk.unit()].push(chunk.bases());
            }
        }
        let factories: Vec<Vec<_>> = contigs
            .iter()
            .map(|&e| (0..e).map(|_| Factory::new()).collect())
            .collect();
        let buffers: Vec<Vec<_>> = contigs
            .iter()
            .map(|&e| (0..e).map(|_| vec![]).collect())
            .collect();
        Self {
            chunks,
            weights,
            factories,
            buffers,
            k,
        }
    }
    pub fn generate_model(
        &mut self,
        ws: &[Vec<f64>],
        reads: &[ERead],
        cl: usize,
        config: &Config,
    ) -> Vec<Vec<DBGHMM>> {
        for contig in self.weights.iter() {
            assert!(!contig.is_empty());
            for unit in contig.iter() {
                assert!(unit.is_empty(), "{:?}", unit);
            }
        }
        for (read, w) in reads.iter().zip(ws) {
            for chunk in read.seq.iter() {
                self.weights[chunk.contig()][chunk.unit()].push(w[cl]);
            }
        }
        let k = self.k;
        assert_eq!(self.weights.len(), self.factories.len());
        assert_eq!(self.weights.len(), self.chunks.len());
        let res: Vec<Vec<_>> = self
            .chunks
            .iter()
            .zip(self.weights.iter())
            .zip(self.factories.iter_mut())
            .zip(self.buffers.iter_mut())
            .map(|(((chunks, weights), fs), bufs)| {
                assert_eq!(chunks.len(), weights.len());
                assert_eq!(chunks.len(), fs.len());
                chunks
                    .par_iter()
                    .zip(weights.par_iter())
                    .zip(fs.par_iter_mut())
                    .zip(bufs.par_iter_mut())
                    .map(|(((cs, w), f), mut buf)| {
                        let mut m = if PRIOR {
                            f.generate_with_weight_prior(&cs, &w, k, &mut buf)
                        } else {
                            f.generate_with_weight(&cs, &w, k)
                        };
                        if MODEL_CHECK {
                            m.check(cs, config, LK_LIMIT);
                        }
                        // if MODEL_CHECK && m.weight() > 2. && m.check(cs, config, LK_LIMIT) {
                        //     // debug!("Model:{} at {} ", m, idx);
                        // }
                        m
                    })
                    .collect()
            })
            .collect();
        for contig in self.weights.iter_mut() {
            for unit in contig.iter_mut() {
                unit.clear();
            }
        }
        res
    }
    pub fn dump_lks(&self, config: &Config, models: &[Vec<Vec<DBGHMM>>]) {
        let max = |x, y| if x < y { y } else { x };
        let min = |x, y| if x < y { x } else { y };
        for (contig, chunks) in self.chunks.iter().enumerate() {
            chunks.par_iter().enumerate().for_each(|(pos, cs)| {
                for c in cs.iter() {
                    let (max, min) = models
                        .iter()
                        .map(|mss| mss[contig][pos].forward(c, config))
                        .fold((-1000000., 10000.), |(x, y), lk| (max(x, lk), min(y, lk)));
                    debug!("LK\t{}\t{}\t{}", contig, pos, max - min);
                }
            })
        }
    }
    fn update_model(
        &mut self,
        ws: &[Vec<f64>],
        mask: &[bool],
        reads: &[ERead],
        cl: usize,
        models: &mut [Vec<DBGHMM>],
        config: &Config,
    ) {
        for ((read, w), &b) in reads.iter().zip(ws).zip(mask) {
            let w = if b { 0. } else { w[cl] };
            for chunk in read.seq.iter() {
                self.weights[chunk.contig()][chunk.unit()].push(w);
            }
        }
        let k = self.k;
        assert_eq!(self.weights.len(), self.factories.len());
        assert_eq!(self.weights.len(), self.chunks.len());
        assert_eq!(self.weights.len(), models.len());
        self.chunks
            .iter()
            .zip(self.weights.iter())
            .zip(self.factories.iter_mut())
            .zip(models.iter_mut())
            .zip(self.buffers.iter_mut())
            .for_each(|((((chunks, weights), fs), ms), bufs)| {
                assert_eq!(chunks.len(), weights.len());
                assert_eq!(chunks.len(), fs.len());
                assert_eq!(chunks.len(), ms.len());
                chunks
                    .par_iter()
                    .zip(bufs.par_iter_mut())
                    .zip(weights.par_iter())
                    .zip(fs.par_iter_mut())
                    .zip(ms.par_iter_mut())
                    .enumerate()
                    .for_each(|(idx, ((((cs, mut buf), w), f), m))| {
                        *m = if PRIOR {
                            // f.generate_with_weight_prior(cs, &w, k, &mut buf)
                            f.update_with_prior(cs, &w, k, &mut buf, m)
                        } else {
                            f.generate_with_weight(cs, &w, k)
                        };
                        if MODEL_CHECK && m.weight() > 2. && m.check(cs, config, LK_LIMIT) {
                            debug!("Model:{} at {}", m, idx);
                        }
                    });
            });
        for contig in self.weights.iter_mut() {
            for unit in contig.iter_mut() {
                unit.clear();
            }
        }
    }
}

fn create_windows(idx: usize, len: usize) -> Vec<(usize, usize, usize)> {
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
        .collect()
}

fn select_within(
    (contig, start, end): (usize, usize, usize),
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    answer: &[u8],
) -> (Vec<ERead>, Vec<u8>, Vec<Vec<u8>>, Vec<u8>) {
    assert_eq!(data.len(), label.len() + answer.len());
    let (contig, start, end) = (contig as u16, start as u16, end as u16);
    let (mut s_data, mut s_label, mut s_forbid, mut s_answer) = (vec![], vec![], vec![], vec![]);
    debug!("Selecting {}\t{}\t{}...", contig, start, end);
    let border = label.len();
    for i in 0..data.len() {
        let read = &data[i];
        let count = read.include_units(contig, start, end);
        if count > 10 {
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
    let nodes_1 = prev.len();
    let nodes_2 = after.len();
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
    let res = bipartite_matching::maximum_weight_matching(nodes_1, nodes_2, &graph);
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
    k: usize,
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
    let clusterings: Vec<Vec<HashSet<String>>> = windows
        .iter()
        .map(|&region| {
            {
                let (contig, start, end) = region;
                debug!("{}:{}-{}", contig, start, end);
            }
            let (data, label, forbidden, answer) =
                select_within(region, data, label, forbidden, answer);
            debug!("Read:{}", data.len());
            assert_eq!(data.len(), forbidden.len());
            assert_eq!(data.len(), label.len() + answer.len());
            let predictions = {
                let (da, la, fo, an) = (&data, &label, &forbidden, &answer);
                clustering(da, la, fo, k, cluster_num, contigs, an, c)
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
    for idx in 0..clusterings.len() {
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
        debug!("Find cluster");
        for (w, clusters) in clusterings.iter().enumerate() {
            for (j, cluster) in clusters.iter().enumerate() {
                if parent == fu.find(j + w * cluster_num).unwrap() {
                    debug!("{}:{}", w, j);
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
                debug!("Read {} does not belong to any cluster.", read.id());
                0
            }
        })
        .collect()
}

pub fn clustering(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
    c: &Config,
) -> Vec<u8> {
    assert_eq!(forbidden.len(), data.len());
    assert_eq!(label.len() + answer.len(), data.len());
    let weights = soft_clustering(data, label, forbidden, k, cluster_num, contigs, answer, c);
    debug!("WEIGHTS\tPrediction. Dump weights");
    assert_eq!(weights.len(), label.len() + answer.len());
    for (weight, ans) in weights.iter().zip(label.iter().chain(answer.iter())) {
        let weights: String = weight
            .iter()
            .map(|e| format!("{:.3},", e))
            .fold(String::new(), |x, y| x + &y);
        debug!("WEIGHTS\t{}\t{}", weights, ans);
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

pub fn entropy(xs: &[f64]) -> f64 {
    assert!(xs.iter().all(|&x| x <= 1.000_000_1 && 0. <= x), "{:?}", xs);
    xs.iter()
        .map(|&x| if x < 0.0001 { 0. } else { -x * x.ln() })
        .sum::<f64>()
}

fn get_max_coverage(data: &[ERead], contig: &[usize]) -> usize {
    let mut cov: Vec<Vec<usize>> = contig.iter().map(|&len| vec![0; len]).collect();
    for read in data {
        for chunk in read.seq().iter() {
            cov[chunk.contig() as usize][chunk.unit() as usize] += 1;
        }
    }
    cov.into_iter()
        .flat_map(|contig| contig.into_iter().max())
        .max()
        .unwrap()
}

#[allow(dead_code)]
fn get_assignments(wor: &[Vec<f64>]) -> Vec<u8> {
    wor.iter()
        .map(|weight| {
            let (cl, _) =
                weight.iter().enumerate().fold(
                    (0, -1.),
                    |(i, m), (j, &w)| if m < w { (j, w) } else { (i, m) },
                );
            cl as u8
        })
        .collect()
}

/// Predict by EM algorithm. the length of return value is the number of test case.
/// The first `label.len()` elements of the `data` should be already classified somehow and
/// the answers should be stored in `label`.
/// When you know the i-th read should not be in the j-th cluster, please add  `j` into `forbidden[i]`'s vector.
/// `cluster_num`should be the number of the cluster.
/// `contigs` should be a map from the index of contig -> number of units.
pub fn soft_clustering(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
    config: &Config,
) -> Vec<Vec<f64>> {
    assert!(cluster_num > 1);
    debug!("{} reads.", data.len());
    debug!("Labels:{}", label.len());
    let informative = forbidden.iter().filter(|e| !e.is_empty()).count();
    debug!("{} informative forbiddens", informative);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let mut weights_of_reads: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let pick_prob = INIT_PICK_PROB;
    let pick_up_len = INIT_PICK_PROB.recip().floor() as usize * LOOP_NUM;
    let max_coverage = get_max_coverage(data, contigs);
    let beta_step = 1. + (BETA_STEP - 1.) / (max_coverage as f64).log10();
    info!("MAX Coverage:{}, Beta step:{:.4}", max_coverage, beta_step);
    let mut beta = search_initial_beta(data, label, forbidden, k, cluster_num, contigs, config);
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, data, cl, config))
        .collect();
    let mut updates = vec![false; data.len()];
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 23);
    let (border, datasize) = (label.len(), data.len() as f64);
    let soe_thr =
        SOE_PER_DATA_ENTROPY * (datasize - border as f64) * (2f64).ln() / (cluster_num as f64).ln();
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
    // mf.dump_lks(config, &models);
    'outer: loop {
        for s in 0.. {
            assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
            assert_eq!(ws.len(), cluster_num);
            let before_soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            let wor = &mut weights_of_reads;
            minibatch_sgd_by(wor, &mut ws, border, data, &models, &updates, beta, config);
            updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
            models.iter_mut().enumerate().for_each(|(cluster, model)| {
                mf.update_model(&weights_of_reads, &updates, data, cluster, model, config);
            });
            let soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            if soe < soe_thr || (before_soe - soe < soe_thr && beta > 1.) {
                break 'outer;
            } else if pick_up_len < s && before_soe < soe {
                break;
            }
        }
        let wr = &weights_of_reads;
        report(id, wr, border, answer, &ws, &models, data, beta, config);
        beta *= (beta_step).max(1.);
    }
    let wr = &weights_of_reads;
    debug!("Finish");
    report(id, wr, border, answer, &ws, &models, data, beta, config);
    weights_of_reads
}

fn construct_initial_weights(
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    data_size: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let border = label.len();
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let num_of_ball = cluster_num * NUM_OF_BALL;
    let denom = (num_of_ball as f64).recip();
    let gen_dist = |idx| {
        let mut choices = vec![true; cluster_num];
        let forbidden: &Vec<u8> = &forbidden[idx + border];
        forbidden
            .iter()
            .for_each(|&cl| choices[cl as usize] = false);
        let choices: Vec<_> = choices
            .into_iter()
            .enumerate()
            .filter_map(|(idx, b)| if b { Some(idx) } else { None })
            .collect();
        let mut bucket = vec![0; cluster_num];
        (0..num_of_ball).for_each(|_| bucket[*choices.choose(&mut rng).unwrap()] += 1);
        bucket.iter().map(|&e| e as f64 * denom).collect::<Vec<_>>()
    };
    let weights: Vec<Vec<_>> = label
        .iter()
        .map(|&e| {
            let mut ws = vec![0.; cluster_num];
            ws[e as usize] = 1.;
            ws
        })
        .chain((0..data_size - border).map(gen_dist))
        .collect();
    assert_eq!(weights.len(), data_size);
    assert!(weights.iter().all(|e| e.len() == cluster_num));
    assert!(weights
        .iter()
        .all(|ws| (ws.iter().sum::<f64>() - 1.).abs() < 0.001));
    weights
}

// Find initial good parameter for \beta.
// It frist starts from INITIAL_BETA, then
// multiple by FACTOR it while SoE after T times sampling would not decrease.
// Then, it divides the beta by FACTOR until SoE after T times sampling would not decrease.
// In other words, it searches for a "platoe" beta for the given data.
fn search_initial_beta(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
) -> f64 {
    let seed = data.iter().map(|e| e.seq.len()).sum::<usize>() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let wor = &weight_of_read;
    let border = label.len();
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut beta = INIT_BETA / FACTOR;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let mut diff = -1.;
    while diff < 0. {
        beta *= FACTOR;
        let c_soe = soe_after_sampling(beta, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    while diff > 0. {
        beta /= FACTOR;
        let c_soe = soe_after_sampling(beta, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    beta
}

fn soe_after_sampling<R: Rng>(
    beta: f64,
    data: &[ERead],
    wor: &[Vec<f64>],
    border: usize,
    rng: &mut R,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let mut updates = vec![false; data.len()];
    // let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl, c))
        .collect();
    // let is_high_quality = find_high_quality(&mut models, data, c);
    let mut wor: Vec<Vec<f64>> = wor.to_vec();
    let pick_prob = INIT_PICK_PROB;
    for _ in 0..SAMPLING {
        updates_flags(&mut updates, &wor, rng, pick_prob, border);
        models.iter_mut().enumerate().for_each(|(cluster, model)| {
            mf.update_model(&wor, &updates, data, cluster, model, c);
        });
        minibatch_sgd_by(&mut wor, &mut ws, border, data, &models, &updates, beta, c);
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn report(
    id: u64,
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    ws: &[f64],
    _models: &[Vec<Vec<DBGHMM>>],
    _data: &[ERead],
    beta: f64,
    _c: &Config,
) -> f64 {
    let correct = weight_of_read
        .iter()
        .skip(border)
        .zip(answer.iter())
        .filter(|&(weights, &ans)| {
            weights
                .iter()
                .filter(|&&g| g > weights[ans as usize])
                .count()
                == 0
        })
        .count();
    let count: Vec<_> = (0..ws.len())
        .map(|cl| {
            weight_of_read
                .iter()
                .filter(|r| r.iter().all(|&w| w <= r[cl]))
                .count()
        })
        .map(|e| format!("{}", e))
        .collect();
    let count = count.join("\t");
    let acc = correct as f64 / answer.len() as f64;
    let ws: Vec<_> = (0..ws.len())
        .map(|cl| weight_of_read.iter().map(|r| r[cl]).sum::<f64>())
        .collect();
    let pi: Vec<_> = ws.iter().map(|e| format!("{:.2}", *e)).collect();
    let pi = pi.join("\t");
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    // for (cl, mss) in _models.iter().enumerate() {
    //     for ms in mss.iter() {
    //         for (pos, m) in ms.iter().enumerate() {
    //             if m.weight() > 5. {
    //                 debug!("{}\t{}\t{}", cl, pos, m);
    //             }
    //         }
    //     }
    // }
    // for (weight, ans) in weight_of_read.iter().skip(border).zip(answer.iter()) {
    //     let weights: String = weight
    //         .iter()
    //         .map(|e| format!("{:.3},", e))
    //         .fold(String::new(), |x, y| x + &y);
    //     debug!("WEIGHTS\t{}\t{}", weights, ans);
    // }
    info!(
        "Summary\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}\t{:.2}",
        id, soe, pi, count, beta, correct, acc
    );
    0.
}

fn updates_flags<R: Rng>(
    updates: &mut [bool],
    weight_of_read: &[Vec<f64>],
    rng: &mut R,
    pick_prob: f64,
    border: usize,
) {
    updates.iter_mut().for_each(|e| *e = false);
    let total = ((weight_of_read.len() - border) as f64 * pick_prob).floor() as usize;
    let mut current = 0;
    let max = weight_of_read.len();
    while current < total {
        let i = rng.gen_range(border, max);
        current += (!updates[i]) as usize;
        updates[i] = true;
    }
}

fn minibatch_sgd_by(
    weight_of_read: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[ERead],
    models: &[Vec<Vec<DBGHMM>>],
    updates: &[bool],
    beta: f64,
    c: &Config,
) {
    let cluster_num = models.len();
    let datasize = data.len() as f64;
    data.iter()
        .zip(weight_of_read.iter_mut())
        .zip(updates.iter())
        .enumerate()
        .skip(border)
        .filter(|&(_, (_, &b))| b)
        .for_each(|(_idx, ((read, weights), _))| {
            compute_log_probs(&models, &ws, &read, weights, c);
            weights.iter_mut().for_each(|w| *w *= beta);
            let tot = utils::logsumexp(weights);
            weights.iter_mut().for_each(|w| *w = (*w - tot).exp());
        });
    ws.iter_mut().enumerate().for_each(|(cl, w)| {
        *w = weight_of_read.iter().map(|e| e[cl]).sum::<f64>() / datasize;
    });
    assert_eq!(ws.len(), cluster_num);
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn compute_log_probs(
    models: &[Vec<Vec<DBGHMM>>],
    ws: &[f64],
    read: &ERead,
    gammas: &mut Vec<f64>,
    c: &Config,
) {
    assert_eq!(models.len(), ws.len());
    assert_eq!(models.len(), gammas.len());
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    models
        .par_iter()
        .zip(gammas.par_iter_mut())
        .zip(ws.par_iter())
        .for_each(|((ms, g), w)| {
            let (lk, num) = read
                .seq
                .par_iter()
                .filter_map(|u| {
                    let model = &ms[u.contig()][u.unit()];
                    if model.weight() < 3. || model.is_broken() {
                        None
                    } else {
                        let lk = model.forward(u.bases(), c);
                        let off = OFFSET * offset(model.weight(), A_PRIOR, B_PRIOR);
                        Some(lk + off)
                    }
                })
                .fold(|| (0., 0.), |(lk, num), x| (lk + x, num + 1.))
                .reduce(|| (0., 0.), |(lk, num), (x, y)| (lk + x, num + y));
            let len = read.seq.len() as f64;
            *g = if num > 1. {
                lk / num * len + w.ln()
            } else {
                LK_LIMIT + w.ln()
            };
        });
    // let weight_thr: Vec<_> = models
    //     .iter()
    //     .map(|model| {
    //         let w = read
    //             .seq
    //             .iter()
    //             .map(|u| model[u.contig()][u.unit()].weight())
    //             .sum::<f64>();
    //         w / read.seq.len() as f64 / 3.
    //     })
    //     .collect();
    // let e_thr = (models.len() as f64).ln() * 0.999;
    // let weights: Vec<Vec<_>> = read
    //     .seq
    //     .par_iter()
    //     .enumerate()
    //     .filter_map(|(idx, u)| {
    //         let lks: Vec<_> = models
    //             .iter()
    //             .zip(weight_thr.iter())
    //             .map(|(model, &thr)| {
    //                 let model = &model[u.contig()][u.unit()];
    //                 if model.weight() < thr || model.is_broken() {
    //                     None
    //                 } else {
    //                     let lk = model.forward(u.bases(), c);
    //                     let offset = OFFSET * offset(model.weight(), A_PRIOR, B_PRIOR);
    //                     Some(lk + offset)
    //                 }
    //             })
    //             .collect();
    //         if lks.iter().any(|e| e.is_none()) {
    //             Some(lks)
    //         } else {
    //             let entropy = calc_entropy(&lks);
    //             if entropy < e_thr {
    //                 Some(lks)
    //             } else {
    //                 // debug!("Removing {} unit({}/{})", idx, entropy, e_thr);
    //                 None
    //             }
    //         }
    //     })
    //     .collect();
    // let len = read.seq.len() as f64;
    // gammas
    //     .iter_mut()
    //     .zip(ws.iter())
    //     .enumerate()
    //     .for_each(|(idx, (g, w))| {
    //         let (lk, num) = weights
    //             .iter()
    //             .filter_map(|ws| ws[idx].as_ref())
    //             .fold((0., 0.), |(lk, num), x| (lk + x, num + 1.));
    //         *g = if num > 0. { lk * len / num } else { LK_LIMIT } + w.ln();
    //     });
    // fn calc_entropy(lks: &[Option<f64>]) -> f64 {
    //     let lks = lks.iter().filter_map(|e| e.as_ref());
    //     let cmp = |a: &&f64, b: &&f64| a.partial_cmp(&b).unwrap();
    //     let max = lks.clone().max_by(cmp).unwrap();
    //     let sum = max + lks.clone().map(|x| (x - max).exp()).sum::<f64>().ln();
    //     lks.clone()
    //         .map(|x| (x - sum).exp())
    //         .map(|x| if x < 0.0001 { 0. } else { -x * x.ln() })
    //         .sum::<f64>()
    // }
}

fn offset(x: f64, a: f64, b: f64) -> f64 {
    (x * a + b).exp()
}

/// Construct DBGHMMs for the `cl`-th cluster.
pub fn construct_with_weights(
    ds: &[ERead],
    gammas: &[Vec<f64>],
    k: usize,
    len: &[usize],
    cl: usize,
) -> Vec<Vec<DBGHMM>> {
    // Contig -> Unit -> Seqs.
    let mut chunks: Vec<Vec<Vec<&[u8]>>> = len.iter().map(|&e| vec![vec![]; e]).collect();
    let mut weights: Vec<Vec<Vec<f64>>> = len.iter().map(|&e| vec![vec![]; e]).collect();
    for (read, ws) in ds.iter().zip(gammas) {
        for chunk in read.seq.iter() {
            chunks[chunk.contig()][chunk.unit()].push(chunk.bases());
            weights[chunk.contig()][chunk.unit()].push(ws[cl]);
        }
    }
    chunks
        .into_par_iter()
        .zip(weights.into_par_iter())
        .map(|(chunks, weights)| {
            let mut f = Factory::new();
            let mut buf = vec![];
            chunks
                .into_iter()
                .zip(weights.into_iter())
                .map(|(cs, ws)| {
                    if PRIOR {
                        f.generate_with_weight_prior(&cs, &ws, k, &mut buf)
                    } else {
                        f.generate_with_weight(&cs, &ws, k)
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Return likelihood of the assignments.
pub fn likelihood_of_assignments(
    data: &[ERead],
    weights_of_reads: &[Vec<f64>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    config: &Config,
) -> f64 {
    assert_eq!(weights_of_reads.len(), data.len());
    assert!(cluster_num > 1);
    let mut mf = ModelFactory::new(contigs, data, k);
    let models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, data, cl, config))
        .collect();
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| weights_of_reads.iter().map(|ws| ws[cl]).sum::<f64>())
        .map(|w| w / data.len() as f64)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    likelihood_of_models(&models, data, &ws, config)
}

fn likelihood_of_models(
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    ws: &[f64],
    c: &Config,
) -> f64 {
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    data.par_iter()
        .map(|read| {
            let gs: Vec<_> = models
                .iter()
                .zip(ws.iter())
                .map(|(model, w)| {
                    let (lk, num) = read
                        .seq
                        .iter()
                        .filter_map(|u| {
                            let model = &model[u.contig()][u.unit()];
                            if model.is_broken() {
                                return None;
                            }
                            let lk = model.forward(u.bases(), c);
                            let offset = OFFSET * offset(model.weight(), A_PRIOR, B_PRIOR);
                            Some(lk + offset)
                        })
                        .fold((0., 0), |(lk, num), x| (lk + x, num + 1));
                    let len = read.seq.len();
                    if num > 0 {
                        lk * len as f64 / num as f64 + w.ln()
                    } else {
                        debug!("Warning:{}:{}, goodModel:{}", read.id(), len, num);
                        LK_LIMIT + w.ln()
                    }
                })
                .collect();
            utils::logsumexp(&gs)
        })
        .sum::<f64>()
}

/// Return the pair of clusters giving the highest gain with
/// respect to likelihood.
/// (cluster number, cluster number, likelihood gain when merging two clusters)
/// The input sequence should be a "weighted" predictions.
pub fn get_mergable_cluster(
    data: &[ERead],
    gammas: &[Vec<f64>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
) -> (f64, u8, u8) {
    let datasize = data.len() as f64;
    let ws: Vec<f64> = gammas
        .iter()
        .map(|g| g.iter().sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    let models: Vec<Vec<_>> = (0..cluster_num)
        .map(|cl| construct_with_weights(data, gammas, k, contigs, cl))
        .collect();
    let before = likelihood_of_models(&models, data, &ws, c);
    let (mut max, mut cluster_a, mut cluster_b) = (std::f64::MIN, 0, 0);
    assert!(cluster_num > 2);
    for i in 0..cluster_num {
        for j in i + 1..cluster_num {
            let lk = likelihood_by_merging(data, &gammas, i, j, cluster_num, k, contigs, c);
            if max < lk {
                cluster_a = i;
                cluster_b = j;
                max = lk;
            }
        }
    }
    (max - before, cluster_a as u8, cluster_b as u8)
}

pub fn likelihood_by_merging(
    data: &[ERead],
    weights_or_reads: &[Vec<f64>],
    i: usize,
    j: usize,
    cluster_num: usize,
    k: usize,
    contigs: &[usize],
    config: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let wor = merge_cluster(&weights_or_reads, i, j, cluster_num);
    let ws: Vec<f64> = (0..cluster_num - 1)
        .map(|cl| wor.iter().map(|ws| ws[cl]).sum::<f64>())
        .map(|w| w / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert!(ws.len() == cluster_num - 1);
    likelihood_of_assignments(data, &wor, k, cluster_num - 1, contigs, config)
}

fn merge_cluster(weights_of_reads: &[Vec<f64>], i: usize, j: usize, cl: usize) -> Vec<Vec<f64>> {
    assert!(i < j);
    weights_of_reads
        .iter()
        .map(|read_weight| {
            let mut ws = vec![0.; cl - 1];
            for (idx, w) in read_weight.iter().enumerate() {
                match idx {
                    x if x < j => ws[idx] += w,
                    x if x == j => ws[i] += w,
                    _ => ws[idx - 1] += w,
                }
            }
            ws
        })
        .collect()
}
