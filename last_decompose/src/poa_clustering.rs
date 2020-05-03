use super::variant_calling;
use super::{ERead, Read};
// use crate::utils::logsumexp;
use poa_hmm::*;
use rand::distributions::Standard;
use rand::{seq::SliceRandom, thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const BETA_INCREASE: f64 = 0.9;
const BETA_DECREASE: f64 = 1.1;
const SMALL_WEIGHT: f64 = 0.000_000_000_0001;
const REP_NUM: usize = 3;
const GIBBS_PRIOR: f64 = 0.02;
const STABLE_LIMIT: u32 = 9;
const IS_STABLE: u32 = 5;
const VARIANT_FRACTION: f64 = 0.90;

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
            .filter(|u| u.contig() < pos.len() && u.unit() < pos[u.contig()].len())
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

fn get_cluster_fraction(asns: &[u8], flags: &[bool], cluster_num: usize) -> Vec<f64> {
    (0..cluster_num)
        .map(|cl| {
            asns.iter()
                .zip(flags)
                .filter(|&(&asn, &b)| !b && asn == cl as u8)
                .count()
        })
        .map(|count| (count as f64).max(SMALL_WEIGHT) / asns.len() as f64)
        .collect()
}

fn get_models<F, R>(
    data: &[Read],
    assignments: &[u8],
    sampled: &[bool],
    (cluster_num, chain_len): (usize, usize),
    rng: &mut R,
    param: (i32, i32, &F),
    use_position: &[bool],
    pick: f64,
) -> Vec<Vec<POA>>
where
    R: Rng,
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let mut chunks: Vec<_> = vec![vec![vec![]; chain_len]; cluster_num];
    let choises: Vec<u8> = (0..cluster_num).map(|e| e as u8).collect();
    for ((read, &asn), &b) in data.iter().zip(assignments.iter()).zip(sampled) {
        if b {
            continue;
        }
        let chosen = *choises
            .choose_weighted(rng, |&k| if k == asn { 1. + pick } else { pick })
            .unwrap();
        for &(pos, ref unit) in read.iter() {
            if use_position[pos] {
                chunks[chosen as usize][pos].push(unit.as_slice());
            }
        }
    }
    let seeds: Vec<_> = rng.sample_iter(Standard).take(chain_len).collect();
    chunks
        .par_iter()
        .map(|cluster| {
            let poa = POA::default();
            cluster
                .par_iter()
                .zip(seeds.par_iter())
                .zip(use_position.par_iter())
                .map(|((cs, &s), &b)| {
                    if b {
                        poa.clone().update(cs, &vec![1.; cs.len()], param, s)
                    } else {
                        poa.clone()
                    }
                })
                .collect()
        })
        .collect()
}

fn update_assignments<R: Rng>(
    models: &[Vec<POA>],
    assignments: &mut [u8],
    data: &[Read],
    sampled: &[bool],
    _rng: &mut R,
    cluster_num: usize,
    betas: &[Vec<Vec<f64>>],
    config: &Config,
    forbidden: &[Vec<u8>],
) -> u32 {
    // let choises: Vec<_> = (0..cluster_num).map(|e| e as u8).collect();
    let fraction_on_positions: Vec<Vec<f64>> = {
        let cl = models[0].len();
        let mut total_count = vec![0; cl];
        let mut counts = vec![vec![0; cl]; cluster_num];
        for (&asn, read) in assignments.iter().zip(data) {
            for &(pos, _) in read.iter() {
                total_count[pos] += 1;
                counts[asn as usize][pos] += 1;
            }
        }
        counts
            .iter()
            .map(|cs| {
                cs.iter()
                    .zip(total_count.iter())
                    .map(|(&c, &t)| c as f64 / t as f64 + SMALL_WEIGHT)
                    .collect()
            })
            .collect()
    };
    assignments
        .iter_mut()
        .zip(data.iter())
        .zip(forbidden.iter())
        .zip(sampled.iter())
        .filter(|&(_, &b)| b)
        .map(|(((asn, read), f), _)| {
            let likelihoods: Vec<(usize, Vec<_>)> = read
                .par_iter()
                .map(|&(pos, ref u)| {
                    let lks = models
                        .par_iter()
                        .map(|ms| ms[pos].forward(u, &config))
                        .collect();
                    (pos, lks)
                })
                .collect();
            let ws: Vec<_> = fraction_on_positions
                .iter()
                .map(|ws| read.iter().map(|&(pos, _)| ws[pos]).sum::<f64>())
                .map(|ws| ws / read.len() as f64)
                .collect();
            let weights: Vec<_> = (0..cluster_num)
                .into_par_iter()
                .map(|l| {
                    (0..cluster_num)
                        .map(|k| {
                            if k != l {
                                let (i, j) = (l.max(k), l.min(k));
                                let prior = ws[k].ln() - ws[l].ln();
                                let lkdiff = likelihoods
                                    .iter()
                                    .map(|&(pos, ref lks)| betas[i][j][pos] * (lks[k] - lks[l]))
                                    .sum::<f64>();
                                prior + lkdiff
                            } else {
                                0.
                            }
                        })
                        .map(|lkdiff| lkdiff.exp())
                        .sum::<f64>()
                        .recip()
                })
                .collect();
            let chosen = {
                let (mut max, mut argmax) = (0., 0);
                for (cl, &p) in weights.iter().enumerate() {
                    if !f.contains(&(cl as u8)) && max < p {
                        max = p;
                        argmax = cl as u8;
                    }
                }
                argmax
            };
            // let line: Vec<_> = weights.iter().map(|e| format!("{:.2}", e)).collect();
            // debug!("[{}]", line.join(","));
            // let chosen = *choises
            //     .choose_weighted(rng, |&k| {
            //         if f.contains(&k) {
            //             0.
            //         } else {
            //             weights[k as usize] + SMALL_WEIGHT
            //         }
            //     })
            //     .unwrap();
            let is_the_same = if *asn == chosen { 0 } else { 1 };
            *asn = chosen;
            is_the_same
        })
        .sum::<u32>()
}

fn add(mut betas: Vec<Vec<Vec<f64>>>, y: Vec<Vec<Vec<f64>>>, rep: usize) -> Vec<Vec<Vec<f64>>> {
    let rep = rep as f64;
    for (clusters, clusters_y) in betas.iter_mut().zip(y) {
        for (bs, bs_y) in clusters.iter_mut().zip(clusters_y) {
            for (b, y) in bs.iter_mut().zip(bs_y) {
                *b += y * y / rep;
            }
        }
    }
    betas
}

fn get_variants<F, R: Rng>(
    data: &[Read],
    asn: &[u8],
    (cluster_num, chain_len): (usize, usize),
    rng: &mut R,
    config: &Config,
    param: (i32, i32, &F),
) -> (Vec<Vec<Vec<f64>>>, f64)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let falses = vec![false; data.len()];
    let init: Vec<Vec<_>> = (0..cluster_num)
        .map(|i| (0..i).map(|_| vec![0.; chain_len]).collect())
        .collect();
    let usepos = vec![true; chain_len];
    let ws = get_cluster_fraction(asn, &falses, cluster_num);
    let (variants, prev_lks) = (0..REP_NUM)
        .map(|_| {
            let tuple = (cluster_num, chain_len);
            let ms = get_models(&data, asn, &falses, tuple, rng, param, &usepos, 0.);
            variant_calling::variant_calling_all_pairs(&ms, &data, config, &ws)
        })
        .fold((init, vec![]), |(xs, mut ps), (y, q)| {
            ps.push(q);
            (add(xs, y, REP_NUM), ps)
        });
    let prev_lk = crate::utils::logsumexp(&prev_lks) - (REP_NUM as f64).ln();
    (variants, prev_lk)
}

pub fn gibbs_sampling<F>(
    data: &[ERead],
    labels: (&[u8], &[u8]),
    f: &[Vec<u8>],
    cluster_num: usize,
    limit: u64,
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
    coverage: usize,
) -> Vec<Option<u8>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    if cluster_num <= 1 || data.len() <= 2 {
        return vec![None; data.len()];
    }
    assert_eq!(f.len(), data.len());
    let per_cluster_coverage = coverage / cluster_num;
    let pick_prob = if per_cluster_coverage < 40 {
        0.005
    } else if per_cluster_coverage < 100 {
        0.045 / 60. * (per_cluster_coverage - 40) as f64 + 0.005
    } else {
        0.05
    };
    let sum = data
        .iter()
        .flat_map(|e| e.id().bytes())
        .map(|e| e as u64)
        .sum::<u64>()
        + 1;
    let params = (limit, pick_prob, sum);
    gibbs_sampling_inner(data, labels, f, cluster_num, params, config, aln).0
}

fn print_lk_gibbs<F>(
    (cluster_num, chain_len): (usize, usize),
    asns: &[u8],
    data: &[Read],
    (label, answer): (&[u8], &[u8]),
    id: u64,
    name: &str,
    param: (i32, i32, &F),
    config: &Config,
) where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let falses = vec![false; data.len()];
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(id);
    let rng = &mut rng;
    let tuple = (cluster_num, chain_len);
    let (variants, _) = get_variants(&data, asns, tuple, rng, config, param);
    let (variants, pos) = select_variants(variants, chain_len);
    let betas = normalize_weights(&variants, 2.);
    let models = get_models(data, asns, &falses, tuple, rng, param, &pos, 0.);
    let fraction_on_positions: Vec<Vec<f64>> = {
        let cl = models[0].len();
        let mut total_count = vec![0; cl];
        let mut counts = vec![vec![0; cl]; cluster_num];
        for (&asn, read) in asns.iter().zip(data) {
            for &(pos, _) in read.iter() {
                total_count[pos] += 1;
                counts[asn as usize][pos] += 1;
            }
        }
        counts
            .iter()
            .map(|cs| {
                cs.iter()
                    .zip(total_count.iter())
                    .map(|(&c, &t)| c as f64 / t as f64 + SMALL_WEIGHT)
                    .collect()
            })
            .collect()
    };
    for (idx, (read, ans)) in data.iter().skip(label.len()).zip(answer).enumerate() {
        let lks = calc_probs(&models, read, config, &betas, &fraction_on_positions);
        trace!(
            "FEATURE\t{}\t{}\t{}\t{}\t{}",
            name,
            id,
            idx,
            ans,
            lks.join("\t")
        );
    }
}

fn gen_assignment<R: Rng>(not_allowed: &[u8], rng: &mut R, c: usize) -> u8 {
    loop {
        let a = rng.gen_range(0, c) as u8;
        if !not_allowed.contains(&a) {
            break a;
        }
    }
}

fn gibbs_sampling_inner<F>(
    data: &[ERead],
    (label, answer): (&[u8], &[u8]),
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    (limit, pick_prob, seed): (u64, f64, u64),
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
) -> (Vec<Option<u8>>, bool)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let param = (aln.ins, aln.del, &aln.score);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let rng = &mut rng;
    let mut assignments: Vec<_> = (0..data.len())
        .map(|idx| {
            if idx < label.len() {
                label[idx]
            } else {
                gen_assignment(&forbidden[idx], rng, cluster_num)
            }
        })
        .collect();
    let mut beta = (cluster_num as f64).powi(2);
    let mut count = 0;
    let mut predictions = std::collections::VecDeque::new();
    let asn = &mut assignments;
    let start = std::time::Instant::now();
    let mut lk = std::f64::NEG_INFINITY;
    let tuple = (cluster_num, chain_len);
    if log_enabled!(log::Level::Trace) {
        print_lk_gibbs(tuple, asn, &data, (label, answer), id, "B", param, config);
    }
    while count < STABLE_LIMIT {
        let (variants, next_lk) = get_variants(&data, asn, tuple, rng, config, param);
        let (variants, pos) = select_variants(variants, chain_len);
        let betas = normalize_weights(&variants, beta);
        beta *= match lk.partial_cmp(&next_lk) {
            Some(std::cmp::Ordering::Less) => BETA_DECREASE,
            Some(std::cmp::Ordering::Greater) => BETA_INCREASE,
            _ => 1.,
        };
        info!("LK\t{}\t{:.3}\t{:.3}\t{:.1}", count, next_lk, lk, beta);
        lk = next_lk;
        let changed_num = (0..pick_prob.recip().ceil() as usize / 2)
            .map(|_| {
                let s: Vec<bool> = (0..data.len())
                    .map(|i| match i.cmp(&label.len()) {
                        std::cmp::Ordering::Less => false,
                        _ => rng.gen_bool(pick_prob),
                    })
                    .collect();
                let ms = get_models(&data, asn, &s, tuple, rng, param, &pos, GIBBS_PRIOR);
                let f = forbidden;
                update_assignments(&ms, asn, &data, &s, rng, cluster_num, &betas, config, f)
            })
            .sum::<u32>();
        debug!("CHANGENUM\t{}", changed_num);
        if changed_num <= (data.len() as f64 * pick_prob / 3.).max(2.) as u32 {
            count += 1;
        } else {
            count = 0;
        }
        predictions.push_back(asn.clone());
        if predictions.len() as u32 > STABLE_LIMIT {
            predictions.pop_front();
            assert_eq!(predictions.len(), STABLE_LIMIT as usize);
        }
        report_gibbs(asn, answer, label, count, id, cluster_num, pick_prob);
        let elapsed = (std::time::Instant::now() - start).as_secs();
        if elapsed > limit && count < STABLE_LIMIT / 2 && predictions.len() as u32 == STABLE_LIMIT {
            info!("Break by timelimit:{:?}", std::time::Instant::now() - start);
            let result = predictions_into_assignments(predictions, cluster_num, data.len());
            return (result, false);
        }
    }
    if log_enabled!(log::Level::Trace) {
        print_lk_gibbs(tuple, asn, &data, (label, answer), id, "A", param, config);
    }
    let result = predictions_into_assignments(predictions, cluster_num, data.len());
    (result, true)
}

fn predictions_into_assignments(
    predictions: std::collections::VecDeque<Vec<u8>>,
    cluster_num: usize,
    data: usize,
) -> Vec<Option<u8>> {
    let maximum_a_posterior = |xs: Vec<u8>| {
        let mut counts: Vec<u32> = vec![0; cluster_num];
        for x in xs {
            counts[x as usize] += 1;
        }
        let (cluster, count): (usize, u32) = counts
            .into_iter()
            .enumerate()
            .max_by_key(|e| e.1)
            .unwrap_or((0, 0));
        if count > IS_STABLE {
            Some(cluster as u8)
        } else {
            None
        }
    };
    predictions
        .into_iter()
        .fold(vec![vec![]; data], |mut acc, xs| {
            assert_eq!(acc.len(), xs.len());
            for (y, x) in acc.iter_mut().zip(xs) {
                y.push(x)
            }
            acc
        })
        .into_iter()
        .map(maximum_a_posterior)
        .collect()
}

fn report_gibbs(asn: &[u8], _ans: &[u8], _lab: &[u8], lp: u32, id: u64, cl: usize, beta: f64) {
    let line: Vec<_> = (0..cl)
        .map(|c| asn.iter().filter(|&&a| a == c as u8).count())
        .map(|e| format!("{}", e))
        .collect();
    info!("Summary\t{}\t{}\t{:.4}\t{}", id, lp, beta, line.join("\t"));
}

fn calc_probs(
    models: &[Vec<POA>],
    read: &Read,
    c: &Config,
    bs: &[Vec<Vec<f64>>],
    fractions: &[Vec<f64>],
) -> Vec<String> {
    let likelihoods: Vec<(usize, Vec<_>)> = read
        .par_iter()
        .map(|&(pos, ref u)| {
            let lks = models.par_iter().map(|ms| ms[pos].forward(u, c)).collect();
            (pos, lks)
        })
        .collect();
    let ws: Vec<_> = fractions
        .iter()
        .map(|ws| read.iter().map(|&(pos, _)| ws[pos]).sum::<f64>())
        .map(|ws| ws / read.len() as f64)
        .collect();
    let cluster_num = ws.len();
    let weights: Vec<_> = (0..cluster_num)
        .into_par_iter()
        .map(|l| {
            (0..cluster_num)
                .map(|k| {
                    if k != l {
                        let (i, j) = (l.max(k), l.min(k));
                        let prior = ws[k].ln() - ws[l].ln();
                        let lkdiff = likelihoods
                            .iter()
                            .map(|&(pos, ref lks)| bs[i][j][pos] * (lks[k] - lks[l]))
                            .sum::<f64>();
                        prior + lkdiff
                    } else {
                        0.
                    }
                })
                .map(|lkdiff| lkdiff.exp())
                .sum::<f64>()
                .recip()
        })
        .collect();
    let sum = weights.iter().sum::<f64>();
    weights.iter().map(|w| format!("{}", w / sum)).collect()
}

fn select_variants(
    mut variants: Vec<Vec<Vec<f64>>>,
    chain: usize,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>) {
    let mut position = vec![false; chain];
    for bss in variants.iter_mut() {
        for bs in bss.iter_mut() {
            let thr = {
                let mut var = bs.clone();
                var.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let pos = (var.len() as f64 * VARIANT_FRACTION).floor() as usize;
                var[pos]
            };
            for (idx, b) in bs.iter_mut().enumerate() {
                if *b < thr {
                    *b = 0.;
                } else {
                    position[idx] = true;
                }
            }
        }
    }
    (variants, position)
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

pub fn predict<F>(
    data: &[ERead],
    labels: &[Option<u8>],
    cluster_num: usize,
    c: &Config,
    input: &[ERead],
    seed: u64,
    aln: &AlnParam<F>,
) -> Vec<Option<u8>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let predictions: std::collections::VecDeque<_> = (0..STABLE_LIMIT)
        .map(|st| {
            let seed = seed + st as u64;
            predict_inner(data, labels, cluster_num, c, input, seed, aln)
        })
        .collect();
    predictions_into_assignments(predictions, cluster_num, input.len())
}

fn predict_inner<F>(
    data: &[ERead],
    labels: &[Option<u8>],
    cluster_num: usize,
    c: &Config,
    input: &[ERead],
    seed: u64,
    aln: &AlnParam<F>,
) -> Vec<u8>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let (data, labels): (Vec<ERead>, Vec<_>) = data
        .into_iter()
        .zip(labels)
        .filter_map(|(d, l)| l.map(|l| (d.clone(), l)))
        .unzip();
    let param = (aln.ins, aln.del, &aln.score);
    let (matrix_pos, chain_len) = to_pos(&data);
    let data = serialize(&data, &matrix_pos);
    let input = serialize(input, &matrix_pos);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let rng = &mut rng;
    let beta = 20.;
    let falses = vec![false; data.len()];
    let tuple = (cluster_num, chain_len);
    let (variants, _) = get_variants(&data, &labels, tuple, rng, c, param);
    let (variants, pos) = select_variants(variants, chain_len);
    let betas = normalize_weights(&variants, beta);
    let ws = get_cluster_fraction(&labels, &falses, cluster_num);
    let tuple = (cluster_num, chain_len);
    let models = get_models(&data, &labels, &falses, tuple, rng, param, &pos, 0.);
    let choises: Vec<_> = (0..cluster_num).map(|e| e as u8).collect();
    input
        .iter()
        .filter_map(|read| {
            let weights: Vec<_> = (0..cluster_num)
                .into_par_iter()
                .map(|l| {
                    let likelihoods: Vec<(usize, Vec<_>)> = read
                        .par_iter()
                        .map(|&(pos, ref u)| {
                            let lks = models.par_iter().map(|ms| ms[pos].forward(u, c)).collect();
                            (pos, lks)
                        })
                        .collect();
                    (0..cluster_num)
                        .map(|k| {
                            if k != l {
                                let (i, j) = (l.max(k), l.min(k));
                                let prior = ws[k].ln() - ws[l].ln();
                                let lkdiff = likelihoods
                                    .iter()
                                    .map(|&(pos, ref lks)| betas[i][j][pos] * (lks[k] - lks[l]))
                                    .sum::<f64>();
                                prior + lkdiff
                            } else {
                                0.
                            }
                        })
                        .map(|lkdiff| lkdiff.exp())
                        .sum::<f64>()
                        .recip()
                })
                .collect();
            choises.choose_weighted(rng, |k| weights[*k as usize]).ok()
        })
        .copied()
        .collect()
}
