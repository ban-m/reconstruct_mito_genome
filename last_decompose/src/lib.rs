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
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
const A: f64 = -0.245;
const B: f64 = 3.6;
const INIT_BETA: f64 = 0.02;
const BETA_STEP: f64 = 1.3;
const MAX_BETA: f64 = 1.;
const INIT_PICK_PROB: f64 = 0.02;
const PICK_PROB_STEP: f64 = 1.05;
const MAX_PICK_PROB: f64 = 0.05;
const MINIMUM_PROB: f64 = 0.001;
const LEARNING_RATE: f64 = 0.1;
const MOMENT: f64 = 0.1;
const LOOP_NUM: usize = 4;
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
        for read in data.into_iter() {
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
        for (read, w) in reads.into_iter().zip(ws) {
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
        for ((read, w), &b) in reads.into_iter().zip(ws).zip(mask) {
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
    let weights = soft_clustering(data, label, forbidden, k, cluster_num, contigs, answer);
    // Maybe we should use randomized choose.
    debug!("Prediction. Dump weights");
    for (weight, ans) in weights.iter().zip(answer) {
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
    assert!(xs.iter().all(|&x| x <= 1.0000001 && 0. <= x));
    xs.iter()
        .map(|&x| if x < 0.0001 { 0. } else { -x * x.ln() })
        .sum::<f64>()
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
    let seed = label.iter().sum::<u8>() as u64 + cluster_num as u64 + data.len() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let border = label.len();
    // weight_of_read[i] = "the vector of each cluster for i-th read"
    let mut weight_of_read: Vec<Vec<f64>> =
        construct_initial_weights(label, forbidden, cluster_num, data.len());
    let mut gammas: Vec<Vec<_>> = vec![vec![0.; cluster_num]; data.len()];
    //let mut moments: Vec<Vec<_>> = vec![vec![1. / cluster_num as f64; cluster_num]; data.len()];
    let mut moments: Vec<Vec<_>> = vec![vec![]; data.len()];
    let datasize = data.len() as f64;
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
    let betas = (0..)
        .map(|i| INIT_BETA * BETA_STEP.powi(i))
        .take_while(|&e| e < MAX_BETA)
        .chain(vec![1.]);
    let pick_probs: Vec<_> = (0..)
        .map(|i| INIT_PICK_PROB * PICK_PROB_STEP.powi(i))
        .take_while(|&e| e < MAX_PICK_PROB)
        .map(|pick_prob| (pick_prob, pick_prob.recip().floor() as usize * LOOP_NUM))
        .flat_map(|(pick_prob, num)| vec![pick_prob; num])
        .collect();
    let mut updates = vec![false; data.len()];
    for beta in betas {
        for &pick_prob in &pick_probs {
            let s = std::time::Instant::now();
            updates_flags(&mut updates, &weight_of_read, &mut rng, pick_prob, beta);
            models.iter_mut().enumerate().for_each(|(cluster, model)| {
                mf.update_model(&weight_of_read, &updates, data, cluster, model);
            });
            let s2 = std::time::Instant::now();
            minibatch_sgd_by(
                &mut weight_of_read,
                &mut gammas,
                &mut moments,
                &mut ws,
                border,
                data,
                &models,
                &updates,
                beta,
            );
            let s3 = std::time::Instant::now();
            debug!("Time\t{}\t{}", (s2 - s).as_millis(), (s3 - s2).as_millis());
        }
        report(&weight_of_read, border, answer, &ws, &models, data, beta);
    }
    weight_of_read
}

fn report(
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    ws: &[f64],
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    beta: f64,
) {
    let correct = weight_of_read
        .iter()
        .skip(border)
        .zip(answer.iter())
        .filter(|&(weights, &ans)| weights.iter().all(|&g| g <= weights[ans as usize]))
        .count();
    let acc = correct as f64 / answer.len() as f64;
    let pi: Vec<_> = ws.iter().map(|e| format!("{:.4}", e)).collect();
    let pi = pi.join("\t");
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let lk = likelihood_of_models(&models, data, &ws);
    info!(
        "Summary\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}",
        lk, soe, pi, beta, correct, acc
    );
}

// let res: Vec<(f64, f64)> = data
//     .par_iter()
//     .zip(weight_of_read.par_iter())
//     .filter_map(|(read, ws)| {
//         ws.iter()
//             .enumerate()
//             .filter_map(|(idx, &e)| if e > 0.99 { Some(idx) } else { None })
//             .nth(0)
//             .map(|idx| (read, idx))
//     })
//     .flat_map(|(read, cl)| {
//         read.seq()
//             .iter()
//             .map(|chunk| {
//                 let m = &models[cl][chunk.contig as usize][chunk.unit as usize];
//                 let lk = m.forward(&chunk.bases(), &DEFAULT_CONFIG) + offset(m.weight(), A, B);
//                 let w = m.weight();
//                 (w, lk)
//             })
//             .collect::<Vec<_>>()
//     })
//     .collect();
// for (w, lk) in res {
//     debug!("DepthLK\t:{}\t{}", w, lk);
// }

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
) {
    let cluster_num = models.len();
    let datasize = data.len() as f64;
    let lr: f64 = LEARNING_RATE * INIT_BETA / beta;
    let ws_gradient = data
        .iter()
        .zip(weight_of_read.iter_mut())
        .zip(gammas.iter_mut())
        .zip(moments.iter_mut())
        .zip(updates.iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .map(|((((read, weights), gamma), moment), _)| {
            compute_log_probs(&models, &ws, &read, gamma);
            gamma.iter_mut().for_each(|g| *g = *g * beta);
            let w = utils::logsumexp(&gamma);
            gamma.iter_mut().for_each(|l| *l = (*l - w).exp());
            assert!((1. - gamma.iter().sum::<f64>()).abs() < 0.001);
            if moment.is_empty() {
                *moment = gamma.clone();
            } else {
                moment.iter_mut().zip(gamma.iter()).for_each(|(m, &g)| {
                    let gradient = g - *m;
                    *m += MOMENT * gradient;
                });
            }
            assert!((1. - moment.iter().sum::<f64>()).abs() < 0.001);
            weights.iter_mut().zip(moment.iter()).for_each(|(w, &m)| {
                let gradient = m - *w;
                *w += lr * gradient;
            });
            assert!((1. - weights.iter().sum::<f64>()).abs() < 0.001);
            assert_eq!(gamma.len(), cluster_num);
            // Convert gamma into moment of PI.
            gamma
                .iter_mut()
                .zip(ws.iter())
                .zip(moment.iter())
                .for_each(|((g, &w), &m)| *g = m - w);
            gamma
        })
        .fold(vec![0.; cluster_num], |mut xs, ys| {
            xs.iter_mut()
                .zip(ys.into_iter())
                .for_each(|(x, y)| *x += *y);
            xs
        });
    assert_eq!(ws_gradient.len(), cluster_num);
    assert!(ws_gradient.iter().sum::<f64>().abs() < 0.0001);
    ws.iter_mut().zip(ws_gradient).for_each(|(w, gradient)| {
        let gradient = gradient / datasize;
        *w += gradient * LEARNING_RATE;
    });
    assert_eq!(ws.len(), cluster_num);
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn compute_log_probs(models: &[Vec<Vec<DBGHMM>>], ws: &[f64], read: &ERead, gammas: &mut Vec<f64>) {
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
                .map(|c| {
                    let model = &model[c.contig()][c.unit()];
                    let lk = model.forward(c.bases(), &DEFAULT_CONFIG);
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
    for (read, ws) in ds.into_iter().zip(gammas) {
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
    let cluster_num = models.len();
    let mut gammas: Vec<_> = vec![vec![0.; cluster_num]; data.len()];
    data.par_iter()
        .zip(gammas.par_iter_mut())
        .map(|(read, gs)| {
            compute_log_probs(&models, &ws, &read, gs);
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
