#[macro_use]
extern crate log;
extern crate bio_utils;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate last_tiling;
extern crate nalgebra as na;
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
use dbg_hmm::{Factory, DBGHMM};
pub use eread::*;
pub use find_breakpoint::initial_clusters;
pub use find_breakpoint::Cluster;
use find_breakpoint::ReadClassify;
use last_tiling::repeat::RepeatPairs;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand::{thread_rng, Rng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::{HashMap, HashSet};
pub mod d3_data;
pub mod error_profile;
pub mod poa_clustering;
use dbg_hmm::Config;
pub use poa_clustering::soft_clustering_poa;
pub mod variant_calling;
const MODEL_CHECK: bool = true;
const NUM_OF_BALL: usize = 100;
const WINDOW_SIZE: usize = 300;
const OVERLAP: usize = 100;
const LK_LIMIT: f64 = -140.;
const OFFSET: f64 = 1.;
const A_PRIOR: f64 = -0.3;
const B_PRIOR: f64 = 6.5;
const INIT_BETA: f64 = 0.008;
// This is the search factor for the doubling-step.
const FACTOR: f64 = 1.4;
// Sampling times for Gibbs sampling.
const SAMPLING: usize = 90;
// This is the factor we multiply at each iteration for beta.
// Note that this factor also scaled for the maximum coverage.
const BETA_STEP: f64 = 1.1;
// Maximum beta. Until this beta, we multiply beta_step for the current beta.
// Loop number for Gibbs sampling.
const LOOP_NUM: usize = 15;
// Initial picking probability.
const INIT_PICK_PROB: f64 = 0.05;
const ENTROPY_THR: f64 = 0.05;
const MIN_WEIGHT: f64 = 0.20;
const K: usize = 6;

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
    config: dbg_hmm::Config,
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
            let crs: Vec<_> = initial_clusters
                .iter()
                .filter(|c| c.is_spanned_by(&read))
                .map(|c| c.id as u8)
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
        initial_clusters.len().max(cluster_num - 1) + 1,
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

pub struct ModelFactory<'a> {
    // Contig -> Unit -> Seqs
    chunks: Vec<Vec<&'a [u8]>>,
    weights: Vec<Vec<f64>>,
    factories: Vec<Factory>,
    buffers: Vec<Vec<f64>>,
    k: usize,
}

impl<'a> ModelFactory<'a> {
    pub fn new(chain_len: usize, data: &'a [Read], k: usize) -> Self {
        let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; chain_len];
        let weights: Vec<Vec<f64>> = vec![vec![]; chain_len];
        for read in data.iter() {
            for x in read.iter() {
                chunks[x.0].push(&x.1);
            }
        }
        let factories: Vec<_> = (0..chain_len).map(|_| Factory::new()).collect();
        let buffers: Vec<_> = vec![vec![]; chain_len];
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
        reads: &[Read],
        cl: usize,
        config: &Config,
    ) -> Vec<DBGHMM> {
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for (read, w) in reads.iter().zip(ws) {
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w[cl]);
            }
        }
        let k = self.k;
        assert_eq!(self.weights.len(), self.factories.len());
        assert_eq!(self.weights.len(), self.chunks.len());
        let res: Vec<_> = self
            .chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .zip(self.factories.par_iter_mut())
            .zip(self.buffers.par_iter_mut())
            .map(|(((chunks, ws), f), mut buf)| {
                let mut m = f.generate_with_weight(chunks, ws, k, &mut buf);
                if MODEL_CHECK {
                    m.check(chunks, config, LK_LIMIT);
                }
                m
            })
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        res
    }
    fn update_model(
        &mut self,
        ws: &[Vec<f64>],
        mask: &[bool],
        reads: &[Read],
        cl: usize,
        models: &mut Vec<DBGHMM>,
        config: &Config,
    ) {
        for ((read, w), &b) in reads.iter().zip(ws).zip(mask) {
            // let w = if b { 0. } else { w[cl] };
            let w = if b { 0. } else { w[cl] + MIN_WEIGHT };
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w);
            }
        }
        let k = self.k;
        assert_eq!(self.weights.len(), self.factories.len());
        assert_eq!(self.weights.len(), self.chunks.len());
        assert_eq!(self.weights.len(), models.len());
        self.chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .zip(self.factories.par_iter_mut())
            .zip(models.par_iter_mut())
            .zip(self.buffers.par_iter_mut())
            .for_each(|((((chunks, ws), f), m), mut buf)| {
                *m = f.generate_with_weight(chunks, ws, k, &mut buf);
                if MODEL_CHECK {
                    m.check(chunks, config, LK_LIMIT);
                }
            });
        self.weights.iter_mut().for_each(|ws| ws.clear());
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
    _k: usize,
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
    cluster_num: usize,
    answer: &[u8],
    c: &Config,
) -> Vec<u8> {
    assert_eq!(forbidden.len(), data.len());
    assert_eq!(label.len() + answer.len(), data.len());
    use poa_clustering::DEFAULT_ALN;
    let weights = soft_clustering_poa(data, label, forbidden, cluster_num, answer, c, &DEFAULT_ALN);
    //let weights = soft_clustering(data, label, forbidden, k, cluster_num, contigs, answer, c);
    debug!("WEIGHTS\tPrediction. Dump weights");
    assert_eq!(weights.len(), label.len() + answer.len());
    for ((read, weight), ans) in data
        .iter()
        .zip(weights.iter())
        .zip(label.iter().chain(answer.iter()))
    {
        let weights: String = weight
            .iter()
            .map(|e| format!("{:.1},", e))
            .fold(String::new(), |x, y| x + &y);
        debug!("WEIGHTS\t{}\t{}\t{}", weights, ans, read.id());
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

fn get_max_coverage(data: &[Read], chain_len: usize) -> usize {
    let mut cov = vec![0; chain_len];
    for &(pos, _) in data.iter().flat_map(|r| r.iter()) {
        cov[pos] += 1;
    }
    cov.into_iter().max().unwrap()
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

// Serialize units in read. In other words,
// We serialize the (contig, unit):(usize, usize) pair into position:usize.
// Return value is (map function, maximum position)
pub fn to_pos(reads: &[ERead]) -> (Vec<Vec<usize>>, usize) {
    let max_contig = reads
        .iter()
        .filter_map(|read| read.seq.iter().map(|e| e.contig()).max())
        .max()
        .unwrap_or(0);
    let minmax_units: Vec<_> = (0..=max_contig)
        .map(|c| {
            let iter = reads
                .iter()
                .flat_map(|read| read.seq.iter())
                .filter(|e| e.contig() == c)
                .map(|e| e.unit());
            let max_unit = iter.clone().max()?;
            let min_unit = iter.clone().min()?;
            Some((min_unit, max_unit))
        })
        .collect();
    let mut res: Vec<_> = minmax_units
        .iter()
        .map(|mm| match mm.as_ref() {
            Some(&(_, max)) => vec![0; max + 1],
            None => vec![],
        })
        .collect();
    let mut len = 0;
    for (contig, mm) in minmax_units.into_iter().enumerate() {
        if let Some((min, max)) = mm {
            for i in min..=max {
                res[contig][i] = len + i - min;
            }
            len += max - min + 1;
        }
    }
    (res, len)
}

fn serialize(data: &[ERead], pos: &[Vec<usize>]) -> Vec<Read> {
    fn serialize_read(read: &ERead, pos: &[Vec<usize>]) -> Read {
        read.seq
            .iter()
            .map(|u| (pos[u.contig()][u.unit()], u.bases().to_vec()))
            .collect()
    }
    data.iter().map(|read| serialize_read(read, pos)).collect()
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
    _contigs: &[usize],
    answer: &[u8],
    config: &Config,
) -> Vec<Vec<f64>> {
    assert!(cluster_num > 1);
    debug!("{} reads.", data.len());
    debug!("Labels:{}", label.len());
    let informative = forbidden.iter().filter(|e| !e.is_empty()).count();
    let (matrix_pos, chain_len) = to_pos(data);
    debug!("MaxLength:{}", chain_len);
    let data = serialize(data, &matrix_pos);
    debug!("{} informative forbiddens", informative);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut weights_of_reads =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    let pick_prob = INIT_PICK_PROB;
    let pick_up_len = INIT_PICK_PROB.recip().floor() as usize * LOOP_NUM;
    let max_coverage = get_max_coverage(&data, chain_len);
    let beta_step = 1. + (BETA_STEP - 1.) / (max_coverage as f64).log10();
    info!("MAX Coverage:{}, Beta step:{:.4}", max_coverage, beta_step);
    let mut beta = search_initial_beta(&data, label, forbidden, k, cluster_num, config, chain_len);
    let mut betas: Vec<_> = vec![beta; chain_len];
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data, k);
    let mut models: Vec<Vec<DBGHMM>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, config))
        .collect();
    let mut updates = vec![false; data.len()];
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 23);
    let (border, datasize) = (label.len(), data.len() as f64);
    let max_entropy = (cluster_num as f64).ln() * datasize;
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
    let mut count = 0;
    'outer: loop {
        for _ in 0..pick_up_len {
            let before_soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            let wor = &mut weights_of_reads;
            minibatch_sgd_by(
                wor, &mut ws, border, &data, &models, &updates, &betas, config,
            );
            updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
            models.iter_mut().enumerate().for_each(|(cluster, model)| {
                mf.update_model(&weights_of_reads, &updates, &data, cluster, model, config);
            });
            let soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            let rate = soe / max_entropy;
            let diff = (before_soe - soe) / max_entropy;
            if rate < ENTROPY_THR || (diff < ENTROPY_THR && beta >= 1.0) {
                break 'outer;
            }
        }
        let lk = likelihood_of_models(&models, &data, &ws, config);
        debug!("DUMP\t{}\t{}\t{}", id, count, lk);
        count += 1;
        let wr = &weights_of_reads;
        report(id, wr, border, answer, &ws, &models, &data, beta, config);
        beta = (beta * beta_step).min(1.);
        let weights = variant_calling::variant_call(&models, &data, config);
        let max = weights
            .iter()
            .map(|&b| b.abs())
            .fold(0., |x, y| if x < y { y } else { x });
        let rate = wr.iter().map(|e| entropy(e)).sum::<f64>() / max_entropy;
        betas = weights
            .iter()
            .map(|b| beta * ((1. - rate) * b.abs() / max + rate))
            .collect();
        // for (idx, (w, b)) in betas.iter().zip(weights).enumerate() {
        //     debug!("VP\t{}\t{:.3}\t{:.3}", idx, w, b);
        // }
    }
    let wr = &weights_of_reads;
    report(id, wr, border, answer, &ws, &models, &data, beta, config);
    weights_of_reads
}

pub fn construct_initial_weights(
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
    data: &[Read],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    c: &Config,
    chain_len: usize,
) -> f64 {
    let seed = data.iter().map(|e| e.len()).sum::<usize>() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let weight_of_read = construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    let border = label.len();
    let mut mf = ModelFactory::new(chain_len, data, k);
    let mut beta = INIT_BETA;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let wor = &weight_of_read;
    loop {
        let betas = &vec![beta; chain_len];
        let next = soe_after_sampling(betas, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        if soe > next {
            break;
        } else {
            debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
            beta *= FACTOR;
        }
    }
    beta /= FACTOR;
    loop {
        let betas = &vec![beta; chain_len];
        let next = soe_after_sampling(betas, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        if soe < next {
            break;
        } else {
            debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
            beta /= FACTOR;
        }
    }
    beta
}

fn soe_after_sampling<R: Rng>(
    betas: &[f64],
    data: &[Read],
    wor: &[Vec<f64>],
    border: usize,
    rng: &mut R,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let mut updates = vec![false; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<DBGHMM>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl, c))
        .collect();
    let mut wor = wor.to_vec();
    for _ in 0..SAMPLING {
        updates_flags(&mut updates, &wor, rng, INIT_PICK_PROB, border);
        minibatch_sgd_by(&mut wor, &mut ws, border, data, &models, &updates, betas, c);
        models.iter_mut().enumerate().for_each(|(cluster, model)| {
            mf.update_model(&wor, &updates, data, cluster, model, c);
        });
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn report(
    id: u64,
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    ws: &[f64],
    _models: &[Vec<DBGHMM>],
    _data: &[Read],
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
    data: &[Read],
    models: &[Vec<DBGHMM>],
    updates: &[bool],
    betas: &[f64],
    c: &Config,
) {
    data.par_iter()
        .zip(weight_of_read.par_iter_mut())
        .zip(updates.par_iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .for_each(|((read, weights), _)| {
            compute_log_probs(models, &ws, read, weights, c, betas);
            let tot = utils::logsumexp(weights);
            weights.iter_mut().for_each(|w| *w = (*w - tot).exp());
            //normalize(weights);
        });
    let datasize = data.len() as f64;
    ws.iter_mut().enumerate().for_each(|(cl, w)| {
        *w = weight_of_read.iter().map(|e| e[cl]).sum::<f64>() / datasize;
    });
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn compute_log_probs(
    models: &[Vec<DBGHMM>],
    ws: &[f64],
    read: &Read,
    weights: &mut Vec<f64>,
    c: &Config,
    betas: &[f64],
) {
    assert_eq!(models.len(), ws.len());
    assert_eq!(models.len(), weights.len());
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    let sum = betas.iter().sum::<f64>();
    models
        .par_iter()
        .zip(weights.par_iter_mut())
        .zip(ws.par_iter())
        .for_each(|((ms, g), w)| {
            let lks: Vec<_> = read
                .par_iter()
                .map(|&(pos, ref u)| {
                    let model = &ms[pos];
                    let beta = betas[pos];
                    if model.weight() < 3. || model.is_broken() {
                        None
                    } else {
                        let lk = model.forward(u, c);
                        let off = OFFSET * offset(model.weight(), A_PRIOR, B_PRIOR);
                        Some(beta * (lk + off))
                    }
                })
                .collect();
            if lks.iter().any(|e| e.is_some()) {
                let ave = lks.iter().filter_map(|e| e.as_ref()).sum::<f64>() / sum;
                let replace = |(idx, &x)| match x {
                    Some(x) => x,
                    None => ave * betas[idx],
                };
                *g = lks.iter().enumerate().map(replace).sum::<f64>() + w.ln();
            } else {
                let beta = betas.iter().sum::<f64>() / sum;
                *g = beta * LK_LIMIT * read.len() as f64 + w.ln();
            }
        });
}

fn offset(x: f64, a: f64, b: f64) -> f64 {
    (x * a + b).exp()
}

/// Construct DBGHMMs for the `cl`-th cluster.
pub fn construct_with_weights(
    ds: &[Read],
    gammas: &[Vec<f64>],
    k: usize,
    cl: usize,
    chain_len: usize,
) -> Vec<DBGHMM> {
    let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; chain_len];
    let mut weights: Vec<Vec<f64>> = vec![vec![]; chain_len];
    for (read, ws) in ds.iter().zip(gammas) {
        for &(pos, ref chunk) in read.iter() {
            chunks[pos].push(chunk);
            weights[pos].push(ws[cl]);
        }
    }
    chunks
        .into_par_iter()
        .zip(weights.into_par_iter())
        .map(|(chunks, weights)| {
            let mut f = Factory::new();
            let mut buf = vec![];
            f.generate_with_weight(&chunks, &weights, k, &mut buf)
        })
        .collect()
}

/// Return likelihood of the assignments.
pub fn likelihood_of_assignments(
    data: &[ERead],
    weights_of_reads: &[Vec<f64>],
    k: usize,
    cluster_num: usize,
    config: &Config,
) -> f64 {
    assert_eq!(weights_of_reads.len(), data.len());
    assert!(cluster_num > 1);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let mut mf = ModelFactory::new(chain_len, &data, k);
    let models: Vec<Vec<DBGHMM>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, config))
        .collect();
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| weights_of_reads.iter().map(|ws| ws[cl]).sum::<f64>())
        .map(|w| w / data.len() as f64)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    likelihood_of_models(&models, &data, &ws, config)
}

fn likelihood_of_models(models: &[Vec<DBGHMM>], data: &[Read], ws: &[f64], c: &Config) -> f64 {
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    data.par_iter()
        .map(|read| {
            let gs: Vec<_> = models
                .iter()
                .zip(ws.iter())
                .map(|(model, w)| {
                    let (lk, num) = read
                        .iter()
                        .filter_map(|&(pos, ref u)| {
                            let model = &model[pos];
                            if model.is_broken() {
                                return None;
                            }
                            let lk = model.forward(u, c);
                            let offset = OFFSET * offset(model.weight(), A_PRIOR, B_PRIOR);
                            Some(lk + offset)
                        })
                        .fold((0., 0), |(lk, num), x| (lk + x, num + 1));
                    let len = read.len();
                    if num > 0 {
                        lk * len as f64 / num as f64 + w.ln()
                    } else {
                        debug!("Warning:{}/{}", len, num);
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
    _contigs: &[usize],
    c: &Config,
) -> (f64, u8, u8) {
    let datasize = data.len() as f64;
    let ws: Vec<f64> = gammas
        .iter()
        .map(|g| g.iter().sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    let before = likelihood_of_assignments(&data, gammas, k, cluster_num, c);
    let (mut max, mut cluster_a, mut cluster_b) = (std::f64::MIN, 0, 0);
    assert!(cluster_num > 2);
    for i in 0..cluster_num {
        for j in i + 1..cluster_num {
            let lk = likelihood_by_merging(&data, &gammas, i, j, cluster_num, k, c);
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
    likelihood_of_assignments(data, &wor, k, cluster_num - 1, config)
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
