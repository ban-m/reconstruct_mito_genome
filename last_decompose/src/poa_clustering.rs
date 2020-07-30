use super::variant_calling;
use super::{ERead, Read};
// use crate::utils::logsumexp;
use poa_hmm::*;
use rand::{seq::SliceRandom, thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
// use rayon::prelude::*;
const BETA_INCREASE: f64 = 1.02;
const BETA_DECREASE: f64 = 1.05;
const BETA_MAX: f64 = 0.8;
const SMALL_WEIGHT: f64 = 0.000_000_001;
const REP_NUM: usize = 3;
const GIBBS_PRIOR: f64 = 0.02;
const STABLE_LIMIT: u32 = 8;
const IS_STABLE: u32 = 3;
const VARIANT_FRACTION: f64 = 0.9;

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
    score: score2,
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

#[allow(clippy::too_many_arguments)]
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
        for &(pos, unit) in read.iter() {
            if use_position[pos] {
                chunks[chosen as usize][pos].push(unit);
            }
        }
    }
    chunks
        .iter_mut()
        .for_each(|cluster| cluster.iter_mut().for_each(|cs| cs.shuffle(rng)));
    chunks
        .iter()
        .map(|cluster| {
            cluster
                .iter()
                .zip(use_position.iter())
                .map(|(cs, &b)| {
                    let len = cs.len().min(30);
                    if b {
                        POA::from_slice(&cs[..len], &vec![1.; len], param)
                    } else {
                        POA::default()
                    }
                })
                .collect()
        })
        .collect()
}

fn get_fraction_on_positions(
    assignments: &[u8],
    cluster_num: usize,
    chain_len: usize,
    data: &[Read],
) -> Vec<Vec<f64>> {
    let mut total_count = vec![0; chain_len];
    let mut counts = vec![vec![0; chain_len]; cluster_num];
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
}

fn get_new_assignment(
    read: &Read,
    fractions: &[Vec<f64>],
    f: &[u8],
    models: &[Vec<POA>],
    betas: &[Vec<Vec<f64>>],
    config: &poa_hmm::Config,
    beta: f64,
) -> u8 {
    // let likelihoods: Vec<(usize, Vec<_>)> = read
    //     .par_iter()
    //     .map(|&(pos, ref u)| {
    //         let lks = models
    //             .par_iter()
    //             .map(|ms| ms[pos].forward(u, &config))
    //             .collect();
    //         (pos, lks)
    //     })
    //     .collect();
    let likelihoods: Vec<(usize, Vec<_>)> = read
        .iter()
        .map(|&(pos, ref u)| {
            let lks = models
                .iter()
                .map(|ms| ms[pos].forward(u, &config))
                .collect();
            (pos, lks)
        })
        .collect();
    let ws: Vec<_> = fractions
        .iter()
        .map(|ws| read.iter().map(|&(pos, _)| ws[pos]).sum::<f64>() / read.len() as f64)
        .collect();
    let cluster_num = fractions.len();
    let weights: Vec<_> = (0..cluster_num)
        .map(|l| {
            (0..cluster_num)
                .map(|k| {
                    if k != l {
                        let (i, j) = (l.max(k), l.min(k));
                        let prior = beta * (ws[k].ln() - ws[l].ln());
                        prior
                            + likelihoods
                                .iter()
                                .map(|&(pos, ref lks)| betas[i][j][pos] * (lks[k] - lks[l]))
                                .sum::<f64>()
                    } else {
                        0.
                    }
                })
                .map(|lkdiff| lkdiff.exp())
                .sum::<f64>()
                .recip()
        })
        .collect();
    let (mut max, mut argmax) = (-0.1, 0);
    for (cl, &p) in weights.iter().enumerate() {
        if !f.contains(&(cl as u8)) && max < p {
            max = p;
            argmax = cl as u8;
        }
    }
    argmax
}

#[allow(clippy::too_many_arguments)]
fn update_assignments(
    models: &[Vec<POA>],
    assignments: &mut [u8],
    data: &[Read],
    sampled: &[bool],
    cluster_num: usize,
    betas: &[Vec<Vec<f64>>],
    config: &Config,
    forbidden: &[Vec<u8>],
    beta: f64,
) -> u32 {
    let fractions: Vec<Vec<f64>> =
        get_fraction_on_positions(assignments, cluster_num, models[0].len(), data);
    assignments
        .iter_mut()
        .zip(data.iter())
        .zip(forbidden.iter())
        .zip(sampled.iter())
        .filter(|&(_, &b)| b)
        .map(|(((asn, read), f), _)| {
            let new_asn = get_new_assignment(read, &fractions, f, &models, &betas, &config, beta);
            let is_the_same = if *asn == new_asn { 0 } else { 1 };
            *asn = new_asn;
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

#[allow(clippy::too_many_arguments)]
pub fn gibbs_sampling<F>(
    data: &[Read],
    labels: &[u8],
    answer: Option<&[u8]>,
    f: &[Vec<u8>],
    chain_len: usize,
    cluster_num: usize,
    limit: u64,
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
    coverage: usize,
) -> Vec<u8>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    if cluster_num <= 1 || data.len() <= 2 {
        return vec![0; data.len()];
    }
    assert_eq!(f.len(), data.len());
    let per_cluster_coverage = coverage / cluster_num;
    let pick_prob = if per_cluster_coverage < 40 {
        0.002
    } else if per_cluster_coverage < 100 {
        0.048 / 60. * (per_cluster_coverage - 40) as f64 + 0.002
    } else {
        0.05
    };
    let sum = data.iter().map(|e| e.len()).sum::<usize>() as u64;
    let params = (limit / 2, pick_prob, sum);
    let (res, ok) = gibbs_sampling_inner(
        data,
        labels,
        answer,
        f,
        chain_len,
        cluster_num,
        params,
        config,
        aln,
    );
    if ok {
        res
    } else {
        let params = (limit / 2, 3. * pick_prob, 3 * sum);
        gibbs_sampling_inner(
            data,
            labels,
            answer,
            f,
            chain_len,
            cluster_num,
            params,
            config,
            aln,
        )
        .0
    }
}

#[allow(clippy::too_many_arguments)]
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

fn gen_assignment<R: Rng>(not_allowed: &[u8], _occed: &[u8], rng: &mut R, c: usize) -> u8 {
    loop {
        let a = rng.gen_range(0, c) as u8;
        if !not_allowed.contains(&a) {
            // && !occed.contains(&a) {
            break a;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn gibbs_sampling_inner<F>(
    data: &[Read],
    label: &[u8],
    answer: Option<&[u8]>,
    forbidden: &[Vec<u8>],
    chain_len: usize,
    cluster_num: usize,
    (limit, pick_prob, seed): (u64, f64, u64),
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
) -> (Vec<u8>, bool)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let param = (aln.ins, aln.del, &aln.score);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let rng = &mut rng;
    let occupied = {
        let mut t = label.to_vec();
        t.sort();
        t.dedup();
        t
    };
    let mut assignments: Vec<_> = (0..data.len())
        .map(|idx| {
            if idx < label.len() {
                label[idx]
            } else {
                gen_assignment(&forbidden[idx], &occupied, rng, cluster_num)
            }
        })
        .collect();
    let mut beta = 0.0005;
    let mut count = 0;
    let mut predictions = std::collections::VecDeque::new();
    let asn = &mut assignments;
    let start = std::time::Instant::now();
    let mut lk = std::f64::NEG_INFINITY;
    let tuple = (cluster_num, chain_len);
    if log_enabled!(log::Level::Trace) {
        if let Some(answer) = answer {
            print_lk_gibbs(tuple, asn, &data, (label, answer), id, "B", param, config);
        }
    }
    while count < STABLE_LIMIT {
        let (variants, next_lk) = get_variants(&data, asn, tuple, rng, config, param);
        let (variants, pos) = select_variants(variants, chain_len);
        let betas = normalize_weights(&variants, 2.);
        beta *= match lk.partial_cmp(&next_lk) {
            Some(std::cmp::Ordering::Less) => BETA_DECREASE,
            Some(std::cmp::Ordering::Greater) => BETA_INCREASE,
            _ => 1.,
        };
        beta = beta.min(BETA_MAX);
        lk = next_lk;
        let changed_num = (0..pick_prob.recip().ceil() as usize)
            .map(|_| {
                let s: Vec<bool> = (0..data.len())
                    .map(|i| match i.cmp(&label.len()) {
                        std::cmp::Ordering::Less => false,
                        _ => rng.gen_bool(pick_prob),
                    })
                    .collect();
                let ms = get_models(&data, asn, &s, tuple, rng, param, &pos, GIBBS_PRIOR);
                let f = forbidden;
                update_assignments(&ms, asn, &data, &s, cluster_num, &betas, config, f, beta)
            })
            .sum::<u32>();
        let has_changed = changed_num <= (data.len() as f64 * 0.01).max(2.) as u32;
        count += has_changed as u32;
        count *= has_changed as u32;
        predictions.push_back(asn.clone());
        if predictions.len() as u32 > STABLE_LIMIT {
            predictions.pop_front();
        }
        report_gibbs(asn, changed_num, count, id, cluster_num, pick_prob);
        let elapsed = (std::time::Instant::now() - start).as_secs();
        if elapsed > limit {
            debug!("Break {} elapsed", elapsed);
            let result = predictions.pop_back().unwrap();
            return (result, false);
        }
    }
    if log_enabled!(log::Level::Trace) {
        if let Some(answer) = answer {
            print_lk_gibbs(tuple, asn, &data, (label, answer), id, "A", param, config);
        }
    }
    let result = predictions.pop_back().unwrap();
    (result, true)
}

// pub fn poa_clustering<F>(
//     data: &[Read],
//     chain_len: usize,
//     label: &[u8],
//     ans: Option<&[u8]>,
//     forbidden: &[Vec<u8>],
//     cluster_num: usize,
//     (limit, pick_prob, seed): (u64, f64, u64),
//     config: &poa_hmm::Config,
//     aln: &AlnParam<F>,
// ) -> Vec<u8>
// where
//     F: Fn(u8, u8) -> i32 + std::marker::Sync,
// {
//     let param = (aln.ins, aln.del, &aln.score);
//     let id: u64 = thread_rng().gen::<u64>() % 100_000;
//     let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
//     let rng = &mut rng;
//     let mut assignments: Vec<_> = (0..data.len())
//         .map(|idx| {
//             if idx < label.len() {
//                 label[idx]
//             } else {
//                 loop {
//                     let a = rng.gen_range(0, cluster_num) as u8;
//                     if !forbidden[idx].contains(&a) {
//                         break a;
//                     }
//                 }
//             }
//         })
//         .collect();
//     let mut beta = 0.0005;
//     let mut count = 0;
//     let mut predictions = std::collections::VecDeque::new();
//     let asn = &mut assignments;
//     let start = std::time::Instant::now();
//     let mut lk = std::f64::NEG_INFINITY;
//     let tuple = (cluster_num, chain_len);
//     if let Some(ans) = ans {
//         if log_enabled!(log::Level::Trace) {
//             print_lk_gibbs(tuple, asn, &data, (label, ans), id, "B", param, config);
//         }
//     }
//     while count < STABLE_LIMIT {
//         let (variants, next_lk) = get_variants(data, asn, tuple, rng, config, param);
//         let (variants, pos) = select_variants(variants, chain_len);
//         let betas = normalize_weights(&variants, 2.);
//         beta *= match lk.partial_cmp(&next_lk) {
//             Some(std::cmp::Ordering::Less) => BETA_DECREASE,
//             Some(std::cmp::Ordering::Greater) => BETA_INCREASE,
//             _ => 1.,
//         };
//         beta = beta.min(BETA_MAX);
//         // info!("LK\t{}\t{:.3}\t{:.3}\t{:.1}", count, next_lk, lk, beta);
//         lk = next_lk;
//         let changed_num = (0..pick_prob.recip().ceil() as usize)
//             .map(|_| {
//                 let s: Vec<bool> = (0..data.len())
//                     .map(|i| match i.cmp(&label.len()) {
//                         std::cmp::Ordering::Less => false,
//                         _ => rng.gen_bool(pick_prob),
//                     })
//                     .collect();
//                 let ms = get_models(&data, asn, &s, tuple, rng, param, &pos, GIBBS_PRIOR);
//                 let f = forbidden;
//                 update_assignments(&ms, asn, &data, &s, cluster_num, &betas, config, f, beta)
//             })
//             .sum::<u32>();
//         let has_changed = changed_num <= (data.len() as f64 * 0.01).max(2.) as u32;
//         count += has_changed as u32;
//         count *= has_changed as u32;
//         report_gibbs(asn, changed_num, count, id, cluster_num, pick_prob);
//         predictions.push_back(asn.clone());
//         if predictions.len() as u32 > STABLE_LIMIT {
//             predictions.pop_front();
//         }
//         let elapsed = (std::time::Instant::now() - start).as_secs();
//         if elapsed > limit {
//             debug!("Break {} elapsed", elapsed);
//             return predictions.pop_back().unwrap();
//         }
//     }
//     if let Some(ans) = ans {
//         if log_enabled!(log::Level::Trace) {
//             print_lk_gibbs(tuple, asn, &data, (label, ans), id, "A", param, config);
//         }
//     }
//     predictions.pop_back().unwrap()
// }

// fn predictions_into_assignments(
//     mut predictions: std::collections::VecDeque<Vec<u8>>,
//     _cluster_num: usize,
//     _data: usize,
// ) -> Vec<Option<u8>> {
// let maximum_a_posterior = |xs: Vec<u8>| {
//     let mut counts: Vec<u32> = vec![0; cluster_num];
//     for x in xs {
//         counts[x as usize] += 1;
//     }
//     let (cluster, count): (usize, u32) = counts
//         .into_iter()
//         .enumerate()
//         .max_by_key(|e| e.1)
//         .unwrap_or((0, 0));
//     if count > IS_STABLE / 2 {
//         Some(cluster as u8)
//     } else {
//         None
//     }
// };
// predictions
//     .pop_back()
//     .unwrap()
//     .into_iter()
//     .map(|x| Some(x))
//     .collect()
// let skip_len = predictions.len() - IS_STABLE as usize;
// predictions
//     .into_iter()
//     .skip(skip_len)
//     .fold(vec![vec![]; data], |mut acc, xs| {
//         assert_eq!(acc.len(), xs.len());
//         for (y, x) in acc.iter_mut().zip(xs) {
//             y.push(x)
//         }
//         acc
//     })
//     .into_iter()
//     .map(maximum_a_posterior)
//     .collect()
// }

fn report_gibbs(asn: &[u8], change_num: u32, lp: u32, id: u64, cl: usize, pp: f64) {
    let line = (0..cl)
        .map(|c| bytecount::count(&asn, c as u8))
        .map(|e| format!("{}", e))
        .collect::<Vec<_>>()
        .join("\t");
    // let mut res = vec![vec![0; cl]; cl];
    // for (&asn, &ans) in asn.iter().zip(lab.iter().chain(ans)) {
    //     res[asn as usize][ans as usize] += 1;
    // }
    // let _res: Vec<_> = res
    //     .iter()
    //     .map(|res| {
    //         let res: Vec<_> = res.iter().map(|x| format!("{:.2}", x)).collect();
    //         res.join("\t")
    //     })
    //     .collect();
    let cn = change_num;
    info!("Summary\t{}\t{}\t{:.4}\t{}\t{}", id, lp, pp, cn, line,);
    //info!("Answers\t{}", res.join("\t"));
}

fn calc_probs(
    models: &[Vec<POA>],
    read: &Read,
    c: &Config,
    bs: &[Vec<Vec<f64>>],
    fractions: &[Vec<f64>],
) -> Vec<String> {
    let likelihoods: Vec<(usize, Vec<_>)> = read
        .iter()
        .map(|&(pos, ref u)| {
            let lks = models.iter().map(|ms| ms[pos].forward(u, c)).collect();
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
