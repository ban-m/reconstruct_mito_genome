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
pub mod utils;
use bio_utils::fasta;
use last_tiling::LastTAB;
use log::Level;
// mod assignments;
mod eread;
mod find_breakpoint;
use dbg_hmm::*;
pub use eread::*;
mod digamma;
use digamma::digamma;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
mod error_profile;
// These A and B are to offset the low-coverage region.
const A: f64 = -0.245;
const B: f64 = 3.6;
// This is the initial reverse temprature.
// Note that, from this rev-temprature, we adjust the start beta by doubling-search.
const INIT_BETA: f64 = 0.02;
// This is the search factor for the doubling-step.
const FACTOR: f64 = 1.4;
// Sampling times for variationl bayes.
const SAMPLING_VB: usize = 10;
// Sampling times for Gibbs sampling.
const SAMPING: usize = 500;
// This is the factor we multiply at each iteration for beta.
// Note that this factor also scaled for the maximum coverage.
const BETA_STEP: f64 = 1.2;
// Maximum beta. Until this beta, we multiply beta_step for the current beta.
const MAX_BETA: f64 = 1.;
// This is the parameter for Diriclet prior.
const ALPHA: f64 = 0.01;
// This is the parameter for de Bruijn prior.
const BETA: f64 = 40.0;
// Loop number for Gibbs sampling.
const LOOP_NUM: usize = 3;
// Loop number for variational Bayes.
const LOOP_NUM_VB: usize = 20;
// Initial picking probability.
const INIT_PICK_PROB: f64 = 0.02;
const PICK_PROB_STEP: f64 = 1.05;
const MAX_PICK_PROB: f64 = 0.05;
const MINIMUM_PROB: f64 = 0.001;
const LEARNING_RATE: f64 = 0.2;
const MOMENT: f64 = 0.2;
const SEED: u64 = 100;
const SOE_PER_DATA_ENTROPY: f64 = 0.05;
const SOE_PER_DATA_ENTROPY_VB: f64 = 0.01;
const K: usize = 6;
/// Main method. Decomposing the reads.
/// You should call "merge" method separatly(?) -- should be integrated with this function.
/// TODO: make a procedure to filter out contained reads.
/// TODO: make a procedure to mask critical regions.
/// TODO: introdce repeats when create critical regions.
pub fn decompose(
    read: Vec<fasta::Record>,
    alignments: Vec<LastTAB>,
    contigs: Vec<fasta::Record>,
) -> Vec<(String, u8)> {
    let config = error_profile::summarize_tab(&alignments, &read, &contigs);
    let contigs = last_tiling::Contigs::new(contigs);
    // Alignment informations are completely (losslessly) encoded into reads.
    let encoded_reads = last_tiling::encoding(&read, &contigs, &alignments);
    // We convert these reads into ERead, a lightweight mode.
    let encoded_reads: Vec<_> = encoded_reads.into_iter().map(ERead::new).collect();
    let critical_regions = critical_regions(&encoded_reads, &contigs, &vec![]);
    if log_enabled!(Level::Debug) {
        for c in &critical_regions {
            debug!("{:?}", c);
        }
    }
    let datasize = encoded_reads.len();
    let mut unassigned_reads: Vec<_> = vec![];
    let mut assigned_reads: Vec<_> = vec![];
    let mut labels: Vec<_> = vec![];
    for read in encoded_reads {
        // TODO
        let matched_cluster = critical_regions
            .iter()
            .enumerate()
            .filter(|(_idx, _cr)| false)
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
    let forbidden = {
        let mut forbidden = vec![vec![]; labels.len()];
        forbidden.extend(unassigned_reads.iter().map(|_read| {
            critical_regions
                .iter()
                .enumerate()
                .filter_map(|(_idx, _cr)| {
                    // TODO
                    // if cr.crossed(read) {
                    //     Some(idx as u8)
                    // } else {n
                    None
                    //}
                })
                .collect::<Vec<u8>>()
        }));
        forbidden
    };
    let dataset: Vec<_> = assigned_reads.into_iter().chain(unassigned_reads).collect();
    let contigs: Vec<_> = (0..contigs.get_num_of_contigs())
        .map(|e| contigs.get_last_unit(e as u16).unwrap() as usize + 1)
        .collect();
    let predicts = clustering(
        &dataset,
        &labels,
        &forbidden,
        K,
        critical_regions.len(),
        &contigs,
        &vec![],
        &config,
    );
    dataset
        .into_iter()
        .zip(labels.iter().chain(predicts.iter()))
        .map(|(read, cl)| (read.id().to_string(), *cl))
        .collect()
}

use std::collections::HashMap;
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
                        let max_sim = cluster
                            .iter()
                            .map(|query| similarity[target][query])
                            .max()
                            .unwrap_or(-1);
                        (idx, max_sim)
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

struct ModelFactory<'a> {
    // Contig -> Unit -> Seqs
    chunks: Vec<Vec<Vec<&'a [u8]>>>,
    weights: Vec<Vec<Vec<f64>>>,
    // Contig -> Unit
    factories: Vec<Vec<Factory>>,
    k: usize,
}

impl<'a> ModelFactory<'a> {
    fn new(contigs: &[usize], data: &'a [ERead], k: usize) -> Self {
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
        Self {
            chunks,
            weights,
            factories,
            k,
        }
    }
    fn generate_model(&mut self, ws: &[Vec<f64>], reads: &[ERead], cl: usize) -> Vec<Vec<DBGHMM>> {
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
            .map(|((chunks, weights), fs)| {
                assert_eq!(chunks.len(), weights.len());
                assert_eq!(chunks.len(), fs.len());
                chunks
                    .par_iter()
                    .zip(weights.par_iter())
                    .zip(fs.par_iter_mut())
                    .map(|((cs, w), f)| f.generate_with_weight(&cs, &w, k))
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
    fn update_model(
        &mut self,
        ws: &[Vec<f64>],
        mask: &[bool],
        reads: &[ERead],
        cl: usize,
        models: &mut [Vec<DBGHMM>],
    ) {
        for ((read, w), &b) in reads.iter().zip(ws).zip(mask) {
            let w = (1 - b as u8) as f64 * w[cl];
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
            .for_each(|(((chunks, weights), fs), ms)| {
                assert_eq!(chunks.len(), weights.len());
                assert_eq!(chunks.len(), fs.len());
                assert_eq!(chunks.len(), ms.len());
                chunks
                    .par_iter()
                    .zip(weights.par_iter())
                    .zip(fs.par_iter_mut())
                    .zip(ms.par_iter_mut())
                    .for_each(|(((cs, w), f), m)| *m = f.generate_with_weight(&cs, &w, k))
            });
        for contig in self.weights.iter_mut() {
            for unit in contig.iter_mut() {
                unit.clear();
            }
        }
    }
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
    let border = label.len();
    assert_eq!(forbidden.len(), data.len());
    let weights = soft_clustering(data, label, forbidden, k, cluster_num, contigs, answer, c);
    // let weights = soft_clustering_full(data, label, forbidden, k, cluster_num, contigs, answer, c);
    //let weights = variational_bayes(data, label, forbidden, k, cluster_num, contigs, answer, c);
    // Maybe we should use randomized choose.
    debug!("Prediction. Dump weights");
    for (weight, ans) in weights.iter().zip(answer).skip(border) {
        let weights: String = weight
            .iter()
            .map(|e| format!("{:.3},", e))
            .fold(String::new(), |x, y| x + &y);
        debug!("{}\t{}", weights, ans);
    }
    weights
        .iter()
        .skip(border)
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

fn entropy(xs: &[f64]) -> f64 {
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

pub fn variational_bayes(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
    config: &Config,
) -> Vec<Vec<f64>> {
    let border = label.len();
    let mut weights_of_reads: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let soe_thr =
        SOE_PER_DATA_ENTROPY_VB * (data.len() - label.len()) as f64 * (cluster_num as f64).ln();
    let mut mf = ModelFactory::new(contigs, data, k);
    // Updates Distributions
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, data, cl))
        .collect();
    let wor = &weights_of_reads;
    let betas = get_schedule(data, contigs, wor, label, k, cluster_num, config);
    let mut alphas: Vec<_> = (0..cluster_num)
        .map(|cl| {
            let n_k = weights_of_reads.iter().map(|e| e[cl]).sum::<f64>();
            (n_k + ALPHA - 1.) * betas[0] + 1.
        })
        .collect();
    let mut log_ros = vec![vec![0.; cluster_num]; data.len()];
    debug!("THR:{}", soe_thr);
    let mut assignments: Vec<u8> = get_assignments(&weights_of_reads);
    for beta in betas {
        let mut has_assignment_changed = true;
        let mut count = 0;
        while has_assignment_changed && count < 50 {
            count += 1;
            // Update weight of reads.
            let wor = &mut weights_of_reads;
            let lr = &mut log_ros;
            batch_vb(wor, lr, &alphas, border, data, &models, beta, config);
            // Updates Distributions
            // There is a "implicit" annealing term at mf.genarate_model.
            models = (0..cluster_num)
                .map(|cl| mf.generate_model(&weights_of_reads, data, cl))
                .collect();
            alphas = (0..cluster_num)
                .map(|cl| weights_of_reads.iter().map(|e| e[cl]).sum::<f64>() + ALPHA)
                .map(|alpha| (alpha - 1.) * beta + 1.)
                .collect();
            let wr = &weights_of_reads;
            report(wr, border, answer, &alphas, &models, data, beta, 1., config);
            let next_assignments = get_assignments(&weights_of_reads);
            has_assignment_changed = assignments
                .iter()
                .zip(next_assignments.iter())
                .any(|(a, b)| a != b);
            assignments = next_assignments;
        }
    }
    weights_of_reads
}

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

fn get_schedule(
    data: &[ERead],
    contigs: &[usize],
    wor: &[Vec<f64>],
    label: &[u8],
    k: usize,
    cluster_num: usize,
    config: &Config,
) -> Vec<f64> {
    let max_coverage = get_max_coverage(data, contigs);
    let beta_step = 1. + (BETA_STEP - 1.) * 2. / (max_coverage as f64).log10();
    info!("MAX Coverage:{}, Beta step:{:.4}", max_coverage, beta_step);
    let init_beta = search_initial_beta_vb(wor, data, label, k, cluster_num, contigs, config);
    info!("Initial Beta {:.4} -> {:.4}", INIT_BETA, init_beta);
    (0..)
        .map(|i| init_beta * beta_step.powi(i - 2))
        .take_while(|&e| e < MAX_BETA)
        .chain(vec![1.])
        .collect()
}

fn search_initial_beta_vb(
    wor: &[Vec<f64>],
    data: &[ERead],
    label: &[u8],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
) -> f64 {
    let weight_of_read = wor.to_vec();
    let wor = &weight_of_read;
    let border = label.len();
    let datasize = data.len() as f64;
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut beta = INIT_BETA;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let mut diff = 0.;
    let thr = SOE_PER_DATA_ENTROPY_VB * datasize * (cluster_num as f64).ln();
    debug!("Start serching initial beta...");
    while diff < thr {
        beta *= FACTOR;
        let c_soe = soe_after_sampling_vb(beta, data, wor, border, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    while diff > thr {
        beta /= FACTOR;
        let c_soe = soe_after_sampling_vb(beta, data, wor, border, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    beta
}

fn soe_after_sampling_vb(
    beta: f64,
    data: &[ERead],
    wor: &[Vec<f64>],
    border: usize,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
) -> f64 {
    let alphas: Vec<_> = (0..cluster_num)
        .map(|cl| wor.iter().map(|e| e[cl]).sum::<f64>() + ALPHA)
        .collect();
    let models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl))
        .collect();
    let mut wor: Vec<Vec<f64>> = wor.to_vec();
    let mut log_ros = vec![vec![0.; cluster_num]; data.len()];
    let lr = &mut log_ros;
    for _ in 0..SAMPLING_VB {
        batch_vb(&mut wor, lr, &alphas, border, data, &models, beta, c);
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn batch_vb(
    weights_of_reads: &mut [Vec<f64>],
    log_ros: &mut [Vec<f64>],
    alphas: &[f64],
    border: usize,
    data: &[ERead],
    models: &[Vec<Vec<DBGHMM>>],
    beta: f64,
    c: &Config,
) {
    let alpha_tot = alphas.iter().sum::<f64>();
    data.par_iter()
        .zip(weights_of_reads.par_iter_mut())
        .zip(log_ros.par_iter_mut())
        .skip(border)
        .for_each(|((read, weights), log_ros)| {
            models
                .iter()
                .zip(alphas.iter())
                .zip(log_ros.iter_mut())
                .for_each(|((ms, &w), log_ro)| {
                    let lk = read
                        .seq
                        .iter()
                        .map(|u| {
                            let model = &ms[u.contig()][u.unit()];
                            let lk = model.forward(u.bases(), c);
                            //let prior = c.null_model(u.bases()) + offset(model.weight(), A, B);
                            //let diff = (prior - lk).exp();
                            // This beta is for scaling/prior de Bruijn graph.
                            //let n_k = w - ALPHA;
                            //let weight = prev_beta * n_k / (BETA + n_k);
                            // lk + (weight + (1. - weight) * diff).ln()
                            lk
                        })
                        .sum::<f64>();
                    let prior = digamma(w) - digamma(alpha_tot);
                    *log_ro = (prior + lk) * beta;
                });
            let log_sum_ro = utils::logsumexp(&log_ros);
            weights
                .iter_mut()
                .zip(log_ros)
                .for_each(|(w, r)| *w = (*r - log_sum_ro).exp());
            assert!((1. - weights.iter().sum::<f64>()).abs() < 0.0001);
        });
}

/// Predict by EM algorithm. the length of return value is the number of test case.
/// The first `label.len()` elements of the `data` should be already classified somehow and
/// the answers should be stored in `label`.
/// When you know the i-th read should not be in the j-th cluster, please add  `j` into `forbidden[i]`'s vector.
/// `cluster_num`should be the number of the cluster.
/// `contigs` should be a map from the index of contig -> number of units.
/// ToDO: Mofify to accept configurations.
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
    // weight_of_read[i] = "the vector of each cluster for i-th read"
    let mut weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let soe_thr =
        SOE_PER_DATA_ENTROPY * (data.len() - label.len()) as f64 * (cluster_num as f64).ln();
    let max_coverage = get_max_coverage(data, contigs);
    let beta_step = 1. + (BETA_STEP - 1.) * 2. / (max_coverage as f64).log10();
    info!("MAX Coverage:{}, Beta step:{:.4}", max_coverage, beta_step);
    let init_beta = search_initial_beta(data, label, forbidden, k, cluster_num, contigs, config);
    info!("Initial Beta {:.4} -> {:.4}", INIT_BETA, init_beta);
    let mut betas: Vec<_> = (0..)
        .map(|i| init_beta * beta_step.powi(i))
        .take_while(|&e| e < MAX_BETA)
        .chain(vec![1.])
        .collect();
    let mut pick_probs: Vec<_> = (0..)
        .map(|i| INIT_PICK_PROB * PICK_PROB_STEP.powi(i))
        .take_while(|&e| e < MAX_PICK_PROB)
        .map(|pick_prob| (pick_prob, pick_prob.recip().floor() as usize * LOOP_NUM))
        .flat_map(|(pick_prob, num)| vec![pick_prob; num])
        .collect();
    let mut lr = LEARNING_RATE;
    for s in 0.. {
        let last_sum_of_entropy = from_weight_of_read(
            &mut weight_of_read,
            data,
            label,
            k,
            cluster_num,
            contigs,
            answer,
            config,
            &betas,
            soe_thr,
            &pick_probs,
            lr,
        );
        if last_sum_of_entropy < soe_thr {
            break;
        } else {
            betas = vec![1.];
            lr = 1.;
            pick_probs
                .iter_mut()
                .for_each(|e| *e *= PICK_PROB_STEP.powi(s as i32));
        }
    }
    weight_of_read
}

pub fn soft_clustering_full(
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
    // weight_of_read[i] = "the vector of each cluster for i-th read"
    let mut weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let max_coverage = get_max_coverage(data, contigs);
    let beta_step = 1. + (BETA_STEP - 1.) * 2. / (max_coverage as f64).log10();
    info!("MAX Coverage:{}, Beta step:{:.4}", max_coverage, beta_step);
    let init_beta =
        search_initial_beta_full(data, label, forbidden, k, cluster_num, contigs, config);
    info!("Initial Beta {:.4} -> {:.4}", INIT_BETA, init_beta);
    let betas: Vec<_> = (0..)
        .map(|i| init_beta * beta_step.powi(i))
        .take_while(|&e| e < MAX_BETA)
        .chain(vec![1.])
        .collect();
    let border = label.len();
    let datasize = data.len() as f64;
    let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    let mut moments: Vec<Vec<_>> = vec![vec![]; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weight_of_read.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert_eq!(ws.len(), cluster_num);
    // Cluster -> Contig -> Unit -> DBG/Vec<CUnit>
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weight_of_read, data, cl))
        .collect();
    for beta in betas {
        let lr = LEARNING_RATE / beta * init_beta;
        for _ in 0..20 {
            models = (0..cluster_num)
                .map(|cl| mf.generate_model(&weight_of_read, data, cl))
                .collect();
            batch_sgd(
                &mut weight_of_read,
                &mut gammas,
                &mut moments,
                &mut ws,
                border,
                data,
                &models,
                beta,
                lr,
                config,
            );
        }
        let wr = &weight_of_read;
        report(&wr, border, answer, &ws, &models, data, beta, lr, config);
    }
    weight_of_read
}

fn from_weight_of_read(
    weight_of_read: &mut [Vec<f64>],
    data: &[ERead],
    label: &[u8],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
    config: &Config,
    betas: &[f64],
    soe_thr: f64,
    pick_probs: &[f64],
    lr: f64,
) -> f64 {
    let seed = label.iter().sum::<u8>() as u64 + cluster_num as u64 + data.len() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let border = label.len();
    let datasize = data.len() as f64;
    let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    let mut moments: Vec<Vec<_>> = vec![vec![]; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weight_of_read.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert_eq!(ws.len(), cluster_num);
    // Cluster -> Contig -> Unit -> DBG/Vec<CUnit>
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weight_of_read, data, cl))
        .collect();
    let mut updates = vec![false; data.len()];
    let mut soe = 100000.; // weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    for &beta in betas {
        let mut soe_diff = 1000000.;
        while soe_diff > soe_thr && soe_thr < soe {
            for &pick_prob in pick_probs {
                updates_flags(&mut updates, &weight_of_read, &mut rng, pick_prob, beta);
                models.iter_mut().enumerate().for_each(|(cluster, model)| {
                    mf.update_model(&weight_of_read, &updates, data, cluster, model);
                });
                minibatch_sgd_by(
                    weight_of_read,
                    &mut gammas,
                    &mut moments,
                    &mut ws,
                    border,
                    data,
                    &models,
                    &updates,
                    beta,
                    lr,
                    config,
                );
            }
            let wr = &weight_of_read;
            report(&wr, border, answer, &ws, &models, data, beta, lr, config);
            let c_soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
            soe_diff = soe - c_soe;
            soe = c_soe;
        }
        if soe < soe_thr || soe_diff < -soe_thr {
            info!("Sum of Entropy {:.4} < {:.4}.", soe, soe_thr);
            info!("Or, it increased. {:.4} -> {:.4}.", soe + soe_diff, soe,);
            info!("Early dropping.");
            break;
        }
    }
    soe
}

const NUM_OF_BALL: usize = 10_000;
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
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(SEED);
    let weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let wor = &weight_of_read;
    let border = label.len();
    let datasize = data.len() as f64;
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut beta = INIT_BETA / FACTOR;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let mut diff = 0.;
    let thr = SOE_PER_DATA_ENTROPY * datasize * (cluster_num as f64).ln();
    while diff < thr {
        beta *= FACTOR;
        let c_soe = soe_after_sampling(beta, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    while diff > thr {
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
    let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    let mut moments: Vec<Vec<_>> = vec![vec![]; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl))
        .collect();
    let mut wor: Vec<Vec<f64>> = wor.to_vec();
    for _ in 0..SAMPING {
        updates_flags(&mut updates, &wor, rng, INIT_PICK_PROB, beta);
        models.iter_mut().enumerate().for_each(|(cluster, model)| {
            mf.update_model(&wor, &updates, data, cluster, model);
        });
        minibatch_sgd_by(
            &mut wor,
            &mut gammas,
            &mut moments,
            &mut ws,
            border,
            data,
            &models,
            &updates,
            beta,
            LEARNING_RATE,
            c,
        );
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn report(
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    ws: &[f64],
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    beta: f64,
    lr: f64,
    c: &Config,
) {
    let correct = weight_of_read
        .iter()
        .skip(border)
        .zip(answer.iter())
        .filter(|&(weights, &ans)| weights.iter().all(|&g| g <= weights[ans as usize]))
        .count();
    let acc = correct as f64 / answer.len() as f64;
    let pi: Vec<_> = ws.iter().map(|e| format!("{:.2}", *e)).collect();
    let pi = pi.join("\t");
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let lk = likelihood_of_models(&models, data, &ws, c);
    info!(
        "Summary\t{:.3}\t{:.3}\t{}\t{:.3}\t{:.3}\t{}\t{:.2}",
        lk, soe, pi, beta, lr, correct, acc
    );
}

fn updates_flags<R: Rng>(
    updates: &mut [bool],
    weight_of_read: &[Vec<f64>],
    rng: &mut R,
    pick_prob: f64,
    beta: f64,
) {
    let datasize = weight_of_read.len() as f64;
    let sum_of_entropy = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    debug!("SoE\t{:.8}\t{:.4}\t{:.4}", sum_of_entropy, pick_prob, beta);
    let denom = sum_of_entropy + MINIMUM_PROB * datasize;
    loop {
        let num_ok = updates
            .iter_mut()
            .zip(weight_of_read.iter())
            .map(|(b, w)| {
                let frac = (entropy(w) + MINIMUM_PROB) / denom;
                let prob = (datasize * pick_prob * frac).min(1.);
                *b = rng.gen_bool(prob);
                *b
            })
            .fold(false, |x, y| x | y);
        if num_ok {
            break;
        }
    }
}

fn minibatch_sgd_by(
    weight_of_read: &mut [Vec<f64>],
    gammas: &mut [Vec<f64>],
    moments: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[ERead],
    models: &[Vec<Vec<DBGHMM>>],
    updates: &[bool],
    beta: f64,
    lr: f64,
    c: &Config,
) {
    let cluster_num = models.len();
    let datasize = data.len() as f64;
    let ws_gradient = data
        .par_iter()
        .zip(weight_of_read.par_iter_mut())
        .zip(gammas.par_iter_mut())
        .zip(moments.par_iter_mut())
        .zip(updates.par_iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .map(|((((read, weights), gamma), moment), _)| {
            compute_log_probs(&models, &ws, &read, gamma, c);
            gamma.iter_mut().for_each(|g| *g *= beta);
            let w = utils::logsumexp(&gamma);
            gamma.iter_mut().for_each(|l| *l = (*l - w).exp());
            debug_assert!((1. - gamma.iter().sum::<f64>()).abs() < 0.001);
            if moment.is_empty() {
                *moment = gamma.clone();
            } else {
                moment.iter_mut().zip(gamma.iter()).for_each(|(m, &g)| {
                    let gradient = g - *m;
                    *m += MOMENT * gradient;
                });
            }
            debug_assert!((1. - moment.iter().sum::<f64>()).abs() < 0.001);
            weights.iter_mut().zip(moment.iter()).for_each(|(w, &m)| {
                let gradient = m - *w;
                *w += lr * gradient;
            });
            debug_assert!((1. - weights.iter().sum::<f64>()).abs() < 0.001);
            assert_eq!(gamma.len(), cluster_num);
            // Convert gamma into moment of PI.
            gamma
                .iter_mut()
                .zip(ws.iter())
                .zip(moment.iter())
                .for_each(|((g, &w), &m)| *g = m - w);
            gamma
        })
        .fold(
            || vec![0.; cluster_num],
            |mut xs, ys| {
                xs.iter_mut().zip(ys.iter()).for_each(|(x, &y)| *x += y);
                xs
            },
        )
        .reduce(
            || vec![0.; cluster_num],
            |mut xs, ys| {
                xs.iter_mut().zip(ys.iter()).for_each(|(x, &y)| *x += y);
                xs
            },
        );
    assert_eq!(ws_gradient.len(), cluster_num);
    debug_assert!(ws_gradient.iter().sum::<f64>().abs() < 0.0001);
    ws.iter_mut().zip(ws_gradient).for_each(|(w, gradient)| {
        let gradient = gradient / datasize;
        *w += gradient * lr;
    });
    assert_eq!(ws.len(), cluster_num);
    debug_assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn search_initial_beta_full(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    c: &Config,
) -> f64 {
    let weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), data.len() as u64);
    let wor = &weight_of_read;
    let border = label.len();
    let datasize = data.len() as f64;
    let mut mf = ModelFactory::new(contigs, data, k);
    let mut beta = INIT_BETA / FACTOR;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let mut diff = 0.;
    let thr = SOE_PER_DATA_ENTROPY * datasize * (cluster_num as f64).ln();
    while diff < thr {
        beta *= FACTOR;
        let c_soe = soe_after_sampling_full(beta, data, wor, border, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    while diff > thr {
        beta /= FACTOR;
        let c_soe = soe_after_sampling_full(beta, data, wor, border, cluster_num, &mut mf, c);
        diff = soe - c_soe;
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, c_soe, diff);
    }
    beta
}

fn soe_after_sampling_full(
    beta: f64,
    data: &[ERead],
    wor: &[Vec<f64>],
    border: usize,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    let mut moments: Vec<Vec<_>> = vec![vec![]; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<Vec<DBGHMM>>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl))
        .collect();
    let mut wor: Vec<Vec<f64>> = wor.to_vec();
    for _ in 0..10 {
        models.iter_mut().enumerate().for_each(|(cluster, model)| {
            mf.update_model(&wor, &vec![false; data.len()], data, cluster, model);
        });
        batch_sgd(
            &mut wor,
            &mut gammas,
            &mut moments,
            &mut ws,
            border,
            data,
            &models,
            beta,
            LEARNING_RATE,
            c,
        );
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn batch_sgd(
    weight_of_read: &mut [Vec<f64>],
    gammas: &mut [Vec<f64>],
    moments: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[ERead],
    models: &[Vec<Vec<DBGHMM>>],
    beta: f64,
    lr: f64,
    c: &Config,
) {
    let cluster_num = models.len();
    let datasize = data.len() as f64;
    let ws_gradient = data
        .par_iter()
        .zip(weight_of_read.par_iter_mut())
        .zip(gammas.par_iter_mut())
        .zip(moments.par_iter_mut())
        .skip(border)
        .map(|(((read, weights), gamma), moment)| {
            compute_log_probs(&models, &ws, &read, gamma, c);
            gamma.iter_mut().for_each(|g| *g *= beta);
            let w = utils::logsumexp(&gamma);
            gamma.iter_mut().for_each(|l| *l = (*l - w).exp());
            debug_assert!((1. - gamma.iter().sum::<f64>()).abs() < 0.001);
            if moment.is_empty() {
                *moment = gamma.clone();
            } else {
                moment.iter_mut().zip(gamma.iter()).for_each(|(m, &g)| {
                    let gradient = g - *m;
                    *m += MOMENT * gradient;
                });
            }
            debug_assert!((1. - moment.iter().sum::<f64>()).abs() < 0.001);
            weights.iter_mut().zip(moment.iter()).for_each(|(w, &m)| {
                let gradient = m - *w;
                *w += lr * gradient;
            });
            debug_assert!((1. - weights.iter().sum::<f64>()).abs() < 0.001);
            assert_eq!(gamma.len(), cluster_num);
            // Convert gamma into moment of PI.
            gamma
                .iter_mut()
                .zip(ws.iter())
                .zip(moment.iter())
                .for_each(|((g, &w), &m)| *g = m - w);
            gamma
        })
        .fold(
            || vec![0.; cluster_num],
            |mut xs, ys| {
                xs.iter_mut().zip(ys.iter()).for_each(|(x, &y)| *x += y);
                xs
            },
        )
        .reduce(
            || vec![0.; cluster_num],
            |mut xs, ys| {
                xs.iter_mut().zip(ys.iter()).for_each(|(x, &y)| *x += y);
                xs
            },
        );
    assert_eq!(ws_gradient.len(), cluster_num);
    debug_assert!(ws_gradient.iter().sum::<f64>().abs() < 0.0001);
    ws.iter_mut().zip(ws_gradient).for_each(|(w, gradient)| {
        let gradient = gradient / datasize;
        *w += gradient * lr;
    });
    assert_eq!(ws.len(), cluster_num);
    debug_assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
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
    models
        .iter()
        .zip(ws.iter())
        .zip(gammas.iter_mut())
        .for_each(|((model, w), g)| {
            *g = read
                .seq
                .par_iter()
                .map(|u| {
                    let model = &model[u.contig()][u.unit()];
                    let lk = model.forward(u.bases(), c);
                    lk + offset(model.weight(), A, B)
                })
                .sum::<f64>()
                + w.ln()
        });
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
            chunks
                .into_iter()
                .zip(weights.into_iter())
                .map(|(cs, ws)| f.generate_with_weight(&cs, &ws, k))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Return likelihood of the assignments.
pub fn likelihood_of_assignments(
    data: &[ERead],
    assignments: &[u8],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    config: &Config,
) -> f64 {
    assert_eq!(assignments.len(), data.len());
    assert!(cluster_num > 1);
    // gammas[i] = "the vector of each cluster for i-th read"
    let forbid: Vec<_> = (0..data.len()).map(|_| vec![]).collect();
    let s = data.len() as u64;
    let gammas: Vec<Vec<f64>> =
        construct_initial_weights(assignments, &forbid, cluster_num, data.len(), s);
    let datasize = data.len() as f64;
    assert_eq!(data.len(), gammas.len());
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| gammas.iter().map(|gs| gs[cl]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    let models: Vec<Vec<_>> = (0..cluster_num)
        .map(|cl| construct_with_weights(data, &gammas, k, contigs, cl))
        .collect();
    likelihood_of_models(&models, data, &ws, config)
}

fn likelihood_of_models(
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    ws: &[f64],
    c: &Config,
) -> f64 {
    let cluster_num = models.len();
    let mut gammas: Vec<_> = vec![vec![0.; cluster_num]; data.len()];
    let sum = ws.iter().sum::<f64>();
    let ws: Vec<_> = ws.iter().map(|&e| e / sum).collect();
    data.par_iter()
        .zip(gammas.par_iter_mut())
        .map(|(read, gs)| {
            compute_log_probs(&models, &ws, &read, gs, c);
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
    gammas: &[Vec<f64>],
    i: usize,
    j: usize,
    cl: usize,
    k: usize,
    contigs: &[usize],
    config: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let gammas = merge_cluster(&gammas, i, j, cl);
    let ws: Vec<f64> = (0..cl - 1)
        .map(|cl| gammas.iter().map(|gs| gs[cl]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert!(ws.len() == cl - 1);
    let models: Vec<Vec<_>> = (0..cl - 1)
        .map(|cl| construct_with_weights(data, &gammas, k, contigs, cl))
        .collect();
    likelihood_of_models(&models, data, &ws, config)
}

fn merge_cluster(gammas: &[Vec<f64>], i: usize, j: usize, cl: usize) -> Vec<Vec<f64>> {
    // Merge weight of j into weight of i
    gammas
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
