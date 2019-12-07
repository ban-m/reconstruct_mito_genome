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
mod assignments;
mod find_breakpoint;
use last_tiling::UNIT_SIZE;
mod eread;
use dbg_hmm::*;
pub use eread::*;
// use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
const A: f64 = -0.2446704;
const B: f64 = 3.6172581;
// const A: f64 = -0.25;
// const B: f64 = 3.0;
pub fn decompose(
    read: Vec<fasta::Record>,
    alignments: Vec<LastTAB>,
    contigs: Vec<fasta::Record>,
) -> Vec<Vec<fasta::Record>> {
    let contigs = last_tiling::Contigs::new(contigs);
    // Alignment informations are completely (losslessly) encoded into reads.
    let encoded_reads = last_tiling::encoding(&read, &contigs, &alignments);
    let critical_regions = critical_regions(&encoded_reads, &contigs);
    if log_enabled!(Level::Debug) {
        for c in critical_regions {
            debug!("{:?}", c);
        }
        return vec![];
    }
    let assignments: Vec<_> = critical_regions
        .into_iter()
        .map(|cr| assignments::local_decompose(&cr, &encoded_reads, &contigs))
        .collect();
    // We no longer need any annotataion for critical region.
    let merge_order = assignments::enumerate_merge_order(&assignments);
    let mut assignments: Vec<_> = assignments.into_iter().map(|e| Some(e)).collect();
    for (from, to) in merge_order {
        // Merge `from` assignment into `to`.
        let f = match &assignments[from] {
            Some(ref f) => f,
            None => panic!("{} should be some value, but none.", from),
        };
        let t = match &assignments[to] {
            Some(ref t) => t,
            None => panic!("{} should be some value, but none.", from),
        };
        let new = assignments::merge_two_assignments(f, t);
        assignments[to] = Some(new);
        assignments[from] = None;
    }
    let mut assignments: Vec<_> = assignments.into_iter().filter_map(|e| e).collect();
    let assignment = assignments.pop().unwrap();
    assert!(assignments.is_empty());
    let mut result = vec![vec![]; assignment.get_num_of_cluster()];
    let mut rng: StdRng = SeedableRng::seed_from_u64(2444);
    // The order of read is critical, since the assignments are just array of weight.
    for (idx, read) in read.into_iter().enumerate() {
        let c = assignment.assign(idx, &mut rng);
        result[c].push(read);
    }
    result
}

use std::collections::HashMap;
pub fn clustering_via_alignment(
    reads: &[usize],
    label: &[u8],
    forbidden: &[Vec<u8>],
    similarity: &[HashMap<usize, i64>],
    cluster_num: usize,
) -> Vec<u8> {
    use rand::Rng;
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

const NUM_OF_BALL: usize = 100;
fn construct_initial_weights(
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    data_size: usize,
) -> Vec<Vec<f64>> {
    let border = label.len();
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data_size as u64);
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

pub fn clustering(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
    answer: &[u8],
) -> Vec<u8> {
    let border = label.len();
    let gammas = soft_clustering(data, label, forbidden, k, cluster_num, contigs, answer);
    // Maybe we should use randomized choose.
    debug!("Prediction. Dump weights");
    for gamma in &gammas {
        let weights: String = gamma
            .iter()
            .map(|e| format!("{:.3},", e))
            .fold(String::new(), |x, y| x + &y);
        debug!("{}", weights);
    }
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(4 * data.len() as u64);
    let choices: Vec<usize> = (0..cluster_num).collect();
    gammas
        .iter()
        .skip(border)
        .map(|gamma| {
            assert_eq!(gamma.len(), choices.len());
            *choices
                .choose_weighted(&mut rng, |&idx| gamma[idx])
                .unwrap() as u8
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
) -> Vec<Vec<f64>> {
    assert!(cluster_num > 1);
    let border = label.len();
    // gammas[i] = "the vector of each cluster for i-th read"
    let mut gammas: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len());
    let datasize = data.len() as f64;
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| gammas.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert_eq!(ws.len(), cluster_num);
    // Cluster -> Contig -> Unit
    let mut models: Vec<Vec<_>> = (0..cluster_num)
        .map(|cl| construct_with_weights(data, &gammas, k, contigs, cl))
        .collect();
    let mut beta = 0.02;
    let step = 1.2;
    let mut lks: Vec<f64> = vec![];
    while beta < 1. {
        // for (contig, &units) in contigs.iter().enumerate() {
        //     for unit in 0..units {
        //         for cl in 0..cluster_num {
        //             debug!("{}\t{}\t{}\t{}", contig, unit, cl, models[cl][contig][unit]);
        //         }
        //     }
        // }
        debug!("Beta:{:.4}", beta);
        for i in 0..40 {
            for (idx, w) in ws.iter().enumerate() {
                debug!("{}\t{:.4}", idx, w);
            }
            let next_lk = data
                .par_iter()
                .zip(gammas.par_iter_mut())
                .skip(border)
                .map(|(read, gamma)| {
                    let mut log_ms = compute_log_probs(&models, &ws, &read);
                    let lk = utils::logsumexp(&log_ms);
                    log_ms.iter_mut().for_each(|log_m| *log_m = *log_m * beta);
                    let w = utils::logsumexp(&log_ms);
                    assert_eq!(gamma.len(), cluster_num);
                    assert_eq!(log_ms.len(), cluster_num);
                    // gamma
                    //     .iter_mut()
                    //     .zip(log_ms.iter())
                    //     .for_each(|(g, l)| *g = (l - w).exp());
                    gamma
                        .iter_mut()
                        .zip(log_ms.iter())
                        .for_each(|(g, l)| *g = (*g + (l - w).exp()) / 2.);
                    assert!((1. - gamma.iter().sum::<f64>()).abs() < 0.001);
                    lk
                })
                .sum::<f64>();
            ws.iter_mut().enumerate().for_each(|(cluster, w)| {
                *w = gammas.iter().map(|g| g[cluster]).sum::<f64>() / datasize
            });
            models.iter_mut().enumerate().for_each(|(cluster, model)| {
                *model = construct_with_weights(data, &gammas, k, contigs, cluster)
            });
            assert_eq!(gammas.len(), answer.len() + label.len());
            let correct = gammas
                .iter()
                .skip(border)
                .zip(answer.iter())
                .filter(|&(gamma, &ans)| gamma.iter().all(|&g| g <= gamma[ans as usize]))
                .count();
            info!("LK\t{:.4}\t{}", next_lk, i);
            info!("Correct\t{}\t{}", correct, answer.len());
            let min: f64 = lks
                .iter()
                .map(|&e| (e - next_lk).abs())
                .fold(1., |x, y| x.min(y));
            if min < 0.001 {
                break;
            } else {
                lks.push(next_lk);
            }
        }
        beta *= step;
    }
    let logms: Vec<_> = data
        .par_iter()
        .map(|read| {
            let xs = compute_log_probs(&models, &ws, read);
            utils::logsumexp(&xs)
        })
        .zip(gammas.par_iter())
        .skip(border)
        .zip(answer.par_iter())
        .map(|((lk, gammas), ans)| {
            let (pred, _) = gammas
                .iter()
                .enumerate()
                .fold(
                    (0, std::f64::MIN),
                    |(x, y), (a, &b)| if y < b { (a, b) } else { (x, y) },
                );
            (lk, pred, ans)
        })
        .collect();
    for (lk, pred, ans) in logms {
        debug!("PREDICTION\t{}\t{}\t{}", lk, pred, ans);
    }
    gammas
}

fn compute_log_probs(models: &[Vec<Vec<DBGHMM>>], ws: &[f64], read: &ERead) -> Vec<f64> {
    assert_eq!(models.len(), ws.len());
    models
        .iter()
        .zip(ws.iter())
        .map(|(model, w)| {
            read.seq
                .iter()
                .map(|c| {
                    let model = &model[c.contig()][c.unit()];
                    let lk = model.forward(c.bases(), &DEFAULT_CONFIG);
                    lk + offset(model.weight(), A, B)
                })
                .sum::<f64>()
                + w.ln()
        })
        .collect()
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
    for (read, ws) in ds.into_iter().zip(gammas) {
        for chunk in read.seq.iter() {
            chunks[chunk.contig()][chunk.unit()].push(chunk.bases());
            weights[chunk.contig()][chunk.unit()].push(ws[cl]);
        }
    }
    chunks
        .into_iter()
        .zip(weights.into_iter())
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
) -> f64 {
    assert_eq!(assignments.len(), data.len());
    assert!(cluster_num > 1);
    // gammas[i] = "the vector of each cluster for i-th read"
    let forbid: Vec<_> = (0..data.len()).map(|_| vec![]).collect();
    let gammas: Vec<Vec<f64>> =
        construct_initial_weights(assignments, &forbid, cluster_num, data.len());
    let datasize = data.len() as f64;
    assert_eq!(data.len(), gammas.len());
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| gammas.iter().map(|gs| gs[cl]).sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    let models: Vec<Vec<_>> = (0..cluster_num)
        .map(|cl| construct_with_weights(data, &gammas, k, contigs, cl))
        .collect();
    likelihood_of_models(&models, data, &ws)
}

fn likelihood_of_models(models: &[Vec<Vec<DBGHMM>>], data: &[ERead], ws: &[f64]) -> f64 {
    data.par_iter()
        .map(|read| {
            let log_ms = compute_log_probs(&models, &ws, &read);
            utils::logsumexp(&log_ms)
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
    let before = likelihood_of_models(&models, data, &ws);
    let (mut max, mut cluster_a, mut cluster_b) = (std::f64::MIN, 0, 0);
    assert!(cluster_num > 2);
    for i in 0..cluster_num {
        for j in i + 1..cluster_num {
            let lk = likelihood_by_merging(data, &gammas, i, j, cluster_num, k, contigs);
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
    likelihood_of_models(&models, data, &ws)
}

fn merge_cluster(gammas: &[Vec<f64>], i: usize, j: usize, cl: usize) -> Vec<Vec<f64>> {
    // Merge weight of j into weight of i
    gammas
        .iter()
        .map(|read_weight| {
            let mut ws = vec![0.; cl - 1];
            for (idx, w) in read_weight.iter().enumerate() {
                if idx < j {
                    ws[idx] += w;
                } else if idx == j {
                    ws[i] += w;
                } else {
                    ws[idx - 1] += w;
                }
            }
            ws
        })
        .collect()
}
