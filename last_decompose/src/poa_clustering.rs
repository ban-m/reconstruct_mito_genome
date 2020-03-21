use super::entropy;
//use super::utils;
use super::variant_calling;
use super::{serialize, to_pos};
use super::{ERead, Read};
use poa_hmm::*;
use rand::{thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const ENTROPY_THR: f64 = 0.35;
const PICK_PROB: f64 = 0.02;
const LRATE: f64 = 0.5;
const BETA_INCREASE: f64 = 1.3;
const BETA_DECREASE: f64 = 1.1;
const INIT_BETA: f64 = 0.2;
// const ENTROPY_STEP: f64 = 0.1;
pub struct AlnParam<F>
where
    F: Fn(u8, u8) -> i32,
{
    ins: i32,
    del: i32,
    score: F,
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        3
    } else {
        -4
    }
}

pub const DEFAULT_ALN: AlnParam<fn(u8, u8) -> i32> = AlnParam {
    ins: -6,
    del: -6,
    score,
};

struct ModelFactory<'a> {
    chunks: Vec<Vec<&'a [u8]>>,
    weights: Vec<Vec<f64>>,
}

impl<'a> ModelFactory<'a> {
    fn new(chain_len: usize, data: &'a [Read]) -> Self {
        let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; chain_len];
        let weights: Vec<Vec<f64>> = vec![vec![]; chain_len];
        for read in data.iter() {
            for x in read.iter() {
                chunks[x.0].push(&x.1);
            }
        }
        Self { chunks, weights }
    }
    fn generate_model<F>(
        &mut self,
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        &AlnParam {
            ins,
            del,
            ref score,
        }: &AlnParam<F>,
    ) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for (read, w) in reads.iter().zip(ws) {
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w[cl]);
            }
        }
        let param = (ins, del, score);
        assert_eq!(self.weights.len(), self.chunks.len());
        let res: Vec<_> = self
            .chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .map(|(chunks, ws)| POA::generate(chunks, ws, param))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        res
    }
    fn update_model<F>(
        &mut self,
        mut ms: Vec<POA>,
        updates: &[bool],
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        &AlnParam {
            ins,
            del,
            ref score,
        }: &AlnParam<F>,
    ) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for ((read, w), &b) in reads.iter().zip(ws).zip(updates) {
            let w = if b { 0. } else { w[cl] };
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w);
            }
        }
        assert_eq!(self.weights.len(), self.chunks.len());
        let parameters = (ins, del, score);
        ms = ms
            .into_par_iter()
            .zip(self.chunks.par_iter())
            .zip(self.weights.par_iter())
            .map(|((m, chunks), ws)| m.update(chunks, ws, parameters))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        ms
    }
}

fn convert(c: &dbg_hmm::Config) -> poa_hmm::Config {
    poa_hmm::Config {
        mismatch: c.mismatch,
        base_freq: c.base_freq,
        p_match: c.p_match,
        p_ins: c.p_ins,
        p_del: c.p_del,
        p_extend_ins: c.p_extend_ins,
        p_extend_del: c.p_extend_del,
        p_del_to_ins: c.p_del_to_ins,
    }
}
use rand::seq::SliceRandom;
pub fn soft_clustering_poa<F>(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    answer: &[u8],
    config: &dbg_hmm::Config,
    aln: &AlnParam<F>,
) -> Vec<Vec<f64>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    assert!(cluster_num > 1);
    let config = &convert(config);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut weights_of_reads =
        crate::construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data);
    debug!("Model factory is built");
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
        .collect();
    debug!("Models have been created");
    let (border, datasize) = (label.len(), data.len() as f64);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 23);
    let pick_num = ((datasize * PICK_PROB).floor() as usize).max(1);
    let mut picks: Vec<_> = (0..data.len()).skip(border).collect();
    picks.shuffle(&mut rng);
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let (mut variants, mut prev_lk): (Vec<Vec<_>>, f64) =
        variant_calling::variant_calling_all_pairs(&models, &data, config, &ws);
    {
        let c = cluster_num as f64;
        variants.iter_mut().for_each(|bss| {
            bss.iter_mut()
                .for_each(|bs| bs.iter_mut().for_each(|b| *b = LRATE * *b + c.recip()))
        });
    };
    if log_enabled!(log::Level::Trace) {
        for (idx, (read, ans)) in data.iter().skip(border).zip(answer).enumerate() {
            let lks = calc_lks(&models, &ws, read, config).join("\t");
            debug!("FEATURE\tBEFORE\t{}\t{}\t{}\t{}", id, idx, ans, lks,);
        }
    }
    report(id, &weights_of_reads, border, answer, cluster_num);
    let mut beta = INIT_BETA;
    let mut count = 0;
    for loop_num in 1.. {
        info!(
            "LK\t{}\t{}\t{:.3}\t{}\t{:.3}",
            id, loop_num, prev_lk, pick_num, beta
        );
        let betas = normalize_weights(&variants, beta);
        while !picks.is_empty() {
            let mut updates = vec![false; data.len()];
            for _ in 0..pick_num {
                if let Some(pos) = picks.pop() {
                    updates[pos] = true;
                }
            }
            let w = &mut weights_of_reads;
            update_weights(w, &mut ws, border, &data, &models, &updates, &betas, config);
            let wor = &weights_of_reads;
            models = models
                .into_iter()
                .enumerate()
                .map(|(cl, m)| mf.update_model(m, &updates, wor, &data, cl, aln))
                .collect();
            if log_enabled!(log::Level::Trace) {
                let lk = variant_calling::get_lk(&models, &data, config, &ws);
                debug!("LK\t{}\t{}\t{}", id, count, lk);
                count += 1;
            }
        }
        report(id, &weights_of_reads, border, answer, cluster_num);
        let (weights, lk) = variant_calling::variant_calling_all_pairs(&models, &data, config, &ws);
        variants.iter_mut().zip(weights).for_each(|(bss, w_bss)| {
            bss.iter_mut().zip(w_bss).for_each(|(bs, ws)| {
                bs.iter_mut()
                    .zip(ws)
                    .for_each(|(b, w)| *b = *b * (1. - LRATE) + w * w * LRATE);
            });
        });
        let soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
        if lk <= prev_lk && soe / datasize / (cluster_num as f64).ln() < ENTROPY_THR {
            info!(
                "LK\t{}\t{}\t{:.3}\t{}\t{:.3}",
                id, loop_num, lk, pick_num, beta
            );
            break;
        } else if lk <= prev_lk {
            beta *= BETA_INCREASE;
        } else {
            beta /= BETA_DECREASE;
        }
        prev_lk = lk;
        picks = (0..data.len()).skip(border).collect();
        picks.shuffle(&mut rng);
    }
    if log_enabled!(log::Level::Trace) {
        for (idx, (read, ans)) in data.iter().skip(border).zip(answer).enumerate() {
            let lks = calc_lks(&models, &ws, read, config).join("\t");
            debug!("FEATURE\tAFTER\t{}\t{}\t{}\t{}", id, idx, ans, lks,);
        }
    }
    weights_of_reads
}

// fn find_initial_beta<F>(
//     beta: f64,
//     variants: &[Vec<Vec<f64>>],
//     data: &[Read],
//     picks: &[usize],
//     (weights_of_reads, ws): (&[Vec<f64>], &[f64]),
//     config: &Config,
//     (border, chain_len, cluster_num): (usize, usize, usize),
//     aln: &AlnParam<F>,
// ) -> f64
// where
//     F: Fn(u8, u8) -> i32 + std::marker::Sync,
// {
//     debug!("SEARCH\t{}", beta);
//     let soe = {
//         let mut ws = ws.to_vec();
//         let mut weights_of_reads = weights_of_reads.to_vec();
//         let datasize = data.len() as f64;
//         let mut picks = picks.to_vec();
//         let pick_num = ((datasize * PICK_PROB).floor() as usize).max(1);
//         let mut mf = ModelFactory::new(chain_len, &data);
//         debug!("Model factory is built");
//         let mut models: Vec<Vec<POA>> = (0..cluster_num)
//             .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
//             .collect();
//         let betas = normalize_weights(variants, beta);
//         while !picks.is_empty() {
//             let mut updates = vec![false; data.len()];
//             for _ in 0..pick_num {
//                 if let Some(pos) = picks.pop() {
//                     updates[pos] = true;
//                 }
//             }
//             let w = &mut weights_of_reads;
//             update_weights(w, &mut ws, border, &data, &models, &updates, &betas, config);
//             let wor = &weights_of_reads;
//             models = models
//                 .into_iter()
//                 .enumerate()
//                 .map(|(cl, m)| mf.update_model(m, &updates, wor, &data, cl, aln))
//                 .collect();
//         }
//         weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>() / datasize
//     };
//     if soe / (cluster_num as f64).ln() > ENTROPY_STEP {
//         find_initial_beta(
//             beta / BETA_DECREASE,
//             variants,
//             data,
//             picks,
//             (weights_of_reads, ws),
//             config,
//             (border, chain_len, cluster_num),
//             aln,
//         )
//     } else {
//         beta
//     }
// }

fn calc_lks(models: &[Vec<POA>], ws: &[f64], read: &Read, c: &Config) -> Vec<String> {
    models
        .iter()
        .zip(ws)
        .map(|(ms, w)| {
            w.ln()
                + read
                    .iter()
                    .map(|&(pos, ref unit)| ms[pos].forward(unit, c))
                    .sum::<f64>()
        })
        .map(|lk| format!("{:.1}", lk))
        .collect()
}

fn normalize_weights(variants: &[Vec<Vec<f64>>], beta: f64) -> Vec<Vec<Vec<f64>>> {
    let max = variants
        .iter()
        .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
        .fold(0., |x, &y| if x < y { y } else { x });
    let mut betas = variants.to_vec();
    betas.iter_mut().for_each(|bss| {
        bss.iter_mut().for_each(|betas| {
            betas.iter_mut().for_each(|b| *b *= beta / max);
        })
    });
    betas
}

fn update_weights(
    weight_of_read: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[Read],
    models: &[Vec<POA>],
    updates: &[bool],
    betas: &[Vec<Vec<f64>>],
    c: &Config,
) {
    data.iter()
        .zip(weight_of_read.iter_mut())
        .zip(updates.iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .for_each(|((read, weights), _)| {
            compute_prior_probability(models, &ws, read, weights, c, betas);
        });
    let datasize = data.len() as f64;
    ws.iter_mut().enumerate().for_each(|(cl, w)| {
        *w = weight_of_read.iter().map(|e| e[cl]).sum::<f64>() / datasize;
    });
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn report(
    id: u64,
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    cluster_num: usize,
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
    let count: Vec<_> = (0..cluster_num)
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
    let ws: Vec<_> = (0..cluster_num)
        .map(|cl| weight_of_read.iter().map(|r| r[cl]).sum::<f64>())
        .collect();
    let pi: Vec<_> = ws.iter().map(|e| format!("{:.2}", *e)).collect();
    let pi = pi.join("\t");
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    info!(
        "Summary\t{}\t{:.3}\t{}\t{}\t{}\t{:.2}",
        id, soe, pi, count, correct, acc,
    );
    acc
}

fn compute_prior_probability(
    models: &[Vec<POA>],
    ws: &[f64],
    read: &Read,
    weights: &mut [f64],
    c: &Config,
    betas: &[Vec<Vec<f64>>],
) {
    let positions: Vec<usize> = read.iter().map(|&(pos, _)| pos).collect();
    let lks: Vec<Vec<f64>> = models
        .par_iter()
        .map(|ms| {
            read.par_iter()
                .map(|&(pos, ref unit)| ms[pos].forward(unit, c))
                .collect()
        })
        .collect();
    let cluster_num = models.len();
    weights.iter_mut().enumerate().for_each(|(k, w)| {
        *w = (0..cluster_num)
            .map(|l| {
                if k == l {
                    0.
                } else {
                    let (i, j) = (l.max(k), l.min(k));
                    let prior = ws[l].ln() - ws[k].ln();
                    let units = positions
                        .iter()
                        .zip(&lks[l])
                        .zip(&lks[k])
                        .map(|((&p, l_lk), k_lk)| betas[i][j][p] * (l_lk - k_lk))
                        .sum::<f64>();
                    prior + units
                }
            })
            .map(|lk_diff| lk_diff.exp())
            .sum::<f64>()
            .recip();
    });
    let sum = weights.iter().sum::<f64>();
    weights.iter_mut().for_each(|w| *w /= sum);
}
