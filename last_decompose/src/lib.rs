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
use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

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

const NUM_OF_BALL: usize = 100;
fn construct_initial_weights_exp(
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

// Construct an initial weights.
// Index of read -> Weights for each cluster.
#[allow(dead_code)]
fn construct_initial_weights(label: &[u8], cluster_num: usize, data_size: usize) -> Vec<Vec<f64>> {
    let border = label.len();
    let uniform = Uniform::new(0, cluster_num as u8);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data_size as u64);
    let init: Vec<u8> = uniform
        .sample_iter(&mut rng)
        .take(data_size - border)
        .collect();
    // For init[i] == i, the weight would be cluster_num + 1/2*cluster_num,
    // Otherwize, the weight would be 1/2*cluster_num;
    let denom = ((2 * cluster_num) as f64).recip();
    let weights: Vec<Vec<_>> = label
        .iter()
        .copied()
        .chain(init)
        .enumerate()
        .map(|(idx, i)| {
            if idx < border {
                let mut ws = vec![0.; cluster_num];
                ws[i as usize] = 1.;
                ws
            } else {
                let mut ws = vec![denom; cluster_num];
                ws[i as usize] = (cluster_num as f64 + 1.) * denom;
                ws
            }
        })
        .collect();
    assert_eq!(weights.len(), data_size);
    assert!(weights.iter().all(|e| e.len() == cluster_num));
    assert!(weights
        .iter()
        .all(|ws| (ws.iter().sum::<f64>() - 1.).abs() < 0.001));
    weights
}

/// Predict by EM algorithm. the length of return value is the number of test case.
/// The first `label.len()` elements of the `data` should be already classified somehow and
/// the answers should be stored in `label`.
/// When you know the i-th read should not be in the j-th cluster, please add  `j` into `forbidden[i]`'s vector.
/// `cluster_num`should be the number of the cluster.
/// `contigs` should be a map from the index of contig -> number of units.
pub fn clustering(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    k: usize,
    cluster_num: usize,
    contigs: &[usize],
) -> Vec<u8> {
    assert!(cluster_num > 1);
    let border = label.len();
    // gammas[i] = "the vector of each cluster for i-th read"
    let mut gammas: Vec<Vec<f64>> =
        construct_initial_weights_exp(label, forbidden, cluster_num, data.len());
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
    for (idx, mss) in models.iter().enumerate() {
        for (contig, ms) in mss.iter().enumerate() {
            for (unit, m) in ms.iter().enumerate() {
                debug!("{}\t{}\t{}\t{}", idx, contig, unit, m);
            }
        }
    }
    let mut beta = 0.03;
    let step = 1.2;
    while beta < 1. {
        debug!("Beta:{:.4}", beta);
        let mut lks: Vec<f64> = vec![];
        for i in 0..30 {
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
                    gamma
                        .iter_mut()
                        .zip(log_ms.iter())
                        .for_each(|(g, l)| *g = (l - w).exp());
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
            info!("LK\t{:.4}\t{}", next_lk, i);
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
    // Maybe we should use randomized choose.
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

fn compute_log_probs(models: &[Vec<Vec<DBGHMM>>], ws: &[f64], read: &ERead) -> Vec<f64> {
    models
        .iter()
        .zip(ws.iter())
        .map(|(model, w)| {
            read.seq
                .iter()
                .map(|c| model[c.contig()][c.unit()].forward(c.bases(), &DEFAULT_CONFIG))
                .sum::<f64>()
                + w.ln()
        })
        .collect()
}

/// Construct DBGHMMs for the `cl`-th cluster.
pub fn construct_with_weights(
    ds: &[ERead],
    gammas: &Vec<Vec<f64>>,
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
            chunks
                .into_par_iter()
                .zip(weights.into_par_iter())
                .map(|(cs, ws)| {
                    let mut f = Factory::new();
                    f.generate_with_weight(&cs, &ws, k)
                })
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
    let gammas: Vec<Vec<f64>> = construct_initial_weights(assignments, cluster_num, data.len());
    let datasize = data.len() as f64;
    let ws: Vec<f64> = gammas
        .iter()
        .map(|g| g.iter().sum::<f64>() / datasize)
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
