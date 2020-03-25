//use super::utils;
use super::variant_calling;
use super::{ERead, Read};
use poa_hmm::*;
use rand::{seq::SliceRandom, thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const ENTROPY_THR: f64 = 0.10;
const PICK_PROB: f64 = 0.01;
const LRATE: f64 = 0.5;
const BETA_INCREASE: f64 = 1.1;
const BETA_DECREASE: f64 = 1.3;
const INIT_BETA: f64 = 0.2;
const NUM_OF_BALL: usize = 100;
const SMALL_WEIGHT: f64 = 0.000_000_000_000_1;
const PRIOR_WEIGHT: f64 = 10.;
// use crate::digamma::digamma;
fn entropy(xs: &[f64]) -> f64 {
    assert!(xs.iter().all(|&x| x <= 1.000_000_1 && 0. <= x), "{:?}", xs);
    xs.iter()
        .map(|&x| if x < 0.0001 { 0. } else { -x * x.ln() })
        .sum::<f64>()
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

pub struct AlnParam<F>
where
    F: Fn(u8, u8) -> i32,
{
    ins: i32,
    del: i32,
    score: F,
}

#[allow(dead_code)]
fn score2(x: u8, y: u8) -> i32 {
    if x == y {
        3
    } else {
        -4
    }
}

#[allow(dead_code)]
fn score(x: u8, y: u8) -> i32 {
    if x == y {
        2
    } else {
        -4
    }
}

pub const DEFAULT_ALN: AlnParam<fn(u8, u8) -> i32> = AlnParam {
    ins: -3,
    del: -3,
    // ins:-6,
    // del:-6,
    score,
};

struct ModelFactory<'a> {
    chunks: Vec<Vec<&'a [u8]>>,
    weights: Vec<Vec<f64>>,
    seeds: Vec<u64>,
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
        let seeds = vec![0; chunks.len()];
        Self {
            chunks,
            weights,
            seeds,
        }
    }
    fn generate_model<F>(
        &mut self,
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        param: &AlnParam<F>,
    ) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        let result = self.chunks.iter().map(|_| POA::default()).collect();
        self.update_model(result, &vec![false; ws.len()], ws, reads, cl, param)
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
        for ((read, w), &_) in reads.iter().zip(ws).zip(updates) {
            // let w = if b { 0. } else { w[cl] };
            let w = w[cl];
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
            .zip(self.seeds.par_iter())
            .map(|(((m, chunks), ws), &s)| m.update(chunks, ws, parameters, s))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        ms
    }
    fn generate_prior_model<F>(&mut self, w: f64, reads: &[Read], aln: &AlnParam<F>) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for read in reads.iter() {
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w);
            }
        }
        assert_eq!(self.weights.len(), self.chunks.len());
        let parameters = (aln.ins, aln.del, &aln.score);
        let ms: Vec<_> = self
            .chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .zip(self.seeds.par_iter())
            .map(|((chunks, ws), &s)| POA::default().update(chunks, ws, parameters, s))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        ms
    }
    fn update_seeds<R: Rng>(&mut self, rng: &mut R) {
        use rand::distributions::Standard;
        self.seeds = rng.sample_iter(Standard).take(self.chunks.len()).collect();
    }
}

pub fn soft_clustering_poa<F>(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    answer: &[u8],
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
) -> Vec<Vec<f64>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    assert!(cluster_num > 1);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 99);
    let mut weights_of_reads =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data);
    mf.update_seeds(&mut rng);
    debug!("Model factory is built");
    // let mut models: Vec<Vec<POA>> = (0..cluster_num)
    //     .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
    //     .collect();
    let updates = vec![false; data.len()];
    let prior_weight = 1. / data.len() as f64;
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| {
            let prior: Vec<POA> = mf.generate_prior_model(prior_weight, &data, aln);
            mf.update_model(prior, &updates, &weights_of_reads, &data, cl, aln)
        })
        .collect();
    debug!("Models have been created");
    let (border, datasize) = (label.len(), data.len() as f64);
    let pick_num = ((datasize * PICK_PROB).floor() as usize).max(1);
    let mut picks: Vec<_> = (0..data.len()).skip(border).collect();
    picks.shuffle(&mut rng);
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let (mut variants, mut prev_lk): (Vec<Vec<_>>, f64) = {
        let models: Vec<Vec<POA>> = (0..cluster_num)
            .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
            .collect();
        let (mut variants, lk) =
            variant_calling::variant_calling_all_pairs(&models, &data, config, &ws);
        let c = cluster_num as f64;
        variants.iter_mut().for_each(|bss| {
            bss.iter_mut()
                .for_each(|bs| bs.iter_mut().for_each(|b| *b = LRATE * *b + c.recip()))
        });
        (variants, lk)
    };
    if log_enabled!(log::Level::Debug) {
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
        // while !picks.is_empty() {
        //     let mut updates = vec![false; data.len()];
        let updates = vec![true; data.len()];
        // for _ in 0..pick_num {
        //     if let Some(pos) = picks.pop() {
        //         updates[pos] = true;
        //     }
        // }
        mf.update_seeds(&mut rng);
        // models = models
        //     .into_iter()
        //     .enumerate()
        //     .map(|(cl, m)| {
        //         mf.update_model(m, &updates, &weights_of_reads, &data, cl,aln)
        //     })
        //     .collect();
        models = (0..cluster_num)
            .map(|cl| {
                mf.update_seeds(&mut rng);
                let prior: Vec<POA> = mf.generate_prior_model(prior_weight, &data, aln);
                mf.update_seeds(&mut rng);
                mf.update_model(prior, &updates, &weights_of_reads, &data, cl, aln)
            })
            .collect();
        let w = &mut weights_of_reads;
        update_weights(
            w, &mut ws, border, &data, &models, &updates, &betas, beta, config,
        );
        if log_enabled!(log::Level::Trace) {
            let lk = variant_calling::get_lk(&models, &data, config, &ws);
            trace!("LK\t{}\t{}\t{}", id, count, lk);
            count += 1;
        }
        //}
        report(id, &weights_of_reads, border, answer, cluster_num);
        let (weights, lk) = {
            let models: Vec<Vec<POA>> = (0..cluster_num)
                .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
                .collect();
            variant_calling::variant_calling_all_pairs(&models, &data, config, &ws)
        };
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
            beta *= BETA_DECREASE;
        } else {
            beta /= BETA_INCREASE;
        }
        prev_lk = lk;
        picks = (0..data.len()).skip(border).collect();
        picks.shuffle(&mut rng);
    }
    if log_enabled!(log::Level::Debug) {
        for (idx, (read, ans)) in data.iter().skip(border).zip(answer).enumerate() {
            let lks = calc_lks(&models, &ws, read, config).join("\t");
            debug!("FEATURE\tAFTER\t{}\t{}\t{}\t{}", id, idx, ans, lks,);
        }
    }
    weights_of_reads
}

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
    beta: f64,
    c: &Config,
) {
    data.iter()
        .zip(weight_of_read.iter_mut())
        .zip(updates.iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .for_each(|((read, weights), _)| {
            compute_prior_probability(models, &ws, read, weights, c, betas, beta);
        });
    ws.iter_mut().enumerate().for_each(|(cl, w)| {
        let new_w = weight_of_read.iter().map(|e| e[cl]).sum::<f64>() / data.len() as f64;
        *w = new_w.max(SMALL_WEIGHT);
    });
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.01, "{:?}", ws);
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
    _beta: f64,
) {
    let positions: Vec<usize> = read.iter().map(|&(pos, _)| pos).collect();
    let lks: Vec<Vec<f64>> = models
        .par_iter()
        .map(|ms| {
            read.par_iter()
                .map(|&(pos, ref unit)| {
                    let m = &ms[pos];
                    m.forward(unit, c) - (1. + PRIOR_WEIGHT / m.weight()).ln()
                })
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
    assert!(
        (1. - weights.iter().sum::<f64>()).abs() < 0.001,
        "{},{:?}",
        sum,
        weights
    );
}

/// Return likelihood of the assignments.
pub fn likelihood_of_assignments<F>(
    data: &[ERead],
    weights_of_reads: &[Vec<f64>],
    cluster_num: usize,
    config: &Config,
    aln: &AlnParam<F>,
) -> f64
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    assert_eq!(weights_of_reads.len(), data.len());
    assert!(cluster_num >= 1);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let mut mf = ModelFactory::new(chain_len, &data);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 9999);
    mf.update_seeds(&mut rng);
    let models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
        .collect();
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| weights_of_reads.iter().map(|ws| ws[cl]).sum::<f64>())
        .map(|w| w / data.len() as f64)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    likelihood_of_models(&models, &data, &ws, config)
}

fn likelihood_of_models(models: &[Vec<POA>], data: &[Read], ws: &[f64], c: &Config) -> f64 {
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    data.par_iter()
        .map(|read| {
            let gs: Vec<_> = models
                .iter()
                .zip(ws.iter())
                .map(|(model, w)| {
                    let lk = read
                        .iter()
                        .map(|&(pos, ref u)| model[pos].forward(u, c))
                        .sum::<f64>();
                    lk + w.ln()
                })
                .collect();
            crate::utils::logsumexp(&gs)
        })
        .sum::<f64>()
}

/// Return the pair of clusters giving the highest gain with
/// respect to likelihood.
/// (cluster number, cluster number, likelihood gain when merging two clusters)
/// The input sequence should be a "weighted" predictions.
pub fn get_mergable_cluster<F>(
    data: &[ERead],
    weights_of_reads: &[Vec<f64>],
    cluster_num: usize,
    c: &Config,
    aln: &AlnParam<F>,
) -> (f64, u8, u8)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let datasize = data.len() as f64;
    let ws: Vec<f64> = weights_of_reads
        .iter()
        .map(|g| g.iter().sum::<f64>() / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    let before = likelihood_of_assignments(&data, weights_of_reads, cluster_num, c, aln);
    let (mut max, mut cluster_a, mut cluster_b) = (std::f64::MIN, 0, 0);
    assert!(cluster_num >= 2);
    for i in 0..cluster_num {
        for j in i + 1..cluster_num {
            let lk = likelihood_by_merging(&data, &weights_of_reads, i, j, cluster_num, c, aln);
            if max < lk {
                cluster_a = i;
                cluster_b = j;
                max = lk;
            }
        }
    }
    let reads_num = weights_of_reads
        .iter()
        .map(|ws| ws[cluster_a] + ws[cluster_b])
        .sum::<f64>();
    let gain_per_read = (max - before) / reads_num;
    (gain_per_read, cluster_a as u8, cluster_b as u8)
}

pub fn likelihood_by_merging<F>(
    data: &[ERead],
    weights_or_reads: &[Vec<f64>],
    i: usize,
    j: usize,
    cluster_num: usize,
    config: &Config,
    aln: &AlnParam<F>,
) -> f64
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let datasize = data.len() as f64;
    let wor = merge_cluster(&weights_or_reads, i, j, cluster_num);
    let ws: Vec<f64> = (0..cluster_num - 1)
        .map(|cl| wor.iter().map(|ws| ws[cl]).sum::<f64>())
        .map(|w| w / datasize)
        .collect();
    assert!((ws.iter().sum::<f64>() - 1.).abs() < 0.0001);
    assert!(ws.len() == cluster_num - 1);
    likelihood_of_assignments(data, &wor, cluster_num - 1, config, aln)
}

pub fn merge_cluster(
    weights_of_reads: &[Vec<f64>],
    i: usize,
    j: usize,
    cl: usize,
) -> Vec<Vec<f64>> {
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
