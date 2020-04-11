use super::variant_calling;
use super::{ERead, Read};
use crate::utils::logsumexp;
use poa_hmm::*;
use rand::distributions::Standard;
use rand::{seq::SliceRandom, thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const ENTROPY_THR: f64 = 0.1;
const BETA_INCREASE: f64 = 0.9;
const BETA_DECREASE: f64 = 1.2;
const INIT_BETA: f64 = 0.2;
const NUM_OF_BALL: usize = 100;
const SMALL_WEIGHT: f64 = 0.000_000_000_0001;
const REP_NUM: usize = 4;
const POA_PRIOR: f64 = 0.5;
const PICK_PRIOR: f64 = 0.3;
const GIBBS_PRIOR: f64 = 0.0;
use crate::digamma::digamma;
const STABLE_LIMIT: u32 = 9;
const IS_STABLE: u32 = 5;
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
    fn update_seeds<R: Rng>(&mut self, rng: &mut R) {
        self.seeds = rng.sample_iter(Standard).take(self.chunks.len()).collect();
    }
}

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
    rng: &mut R,
    cluster_num: usize,
    betas: &[Vec<Vec<f64>>],
    config: &Config,
    forbidden: &[Vec<u8>],
) -> u32 {
    let choises: Vec<_> = (0..cluster_num).map(|e| e as u8).collect();
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
    let seeds: Vec<u64> = rng.sample_iter(Standard).take(data.len()).collect();
    assignments
        .iter_mut()
        .zip(data.iter())
        .zip(forbidden.iter())
        .zip(seeds.into_iter())
        .zip(sampled.iter())
        .filter(|&(_, &b)| b)
        .map(|((((asn, read), f), s), _)| {
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
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(s);
            let chosen = *choises
                .choose_weighted(&mut rng, |k| {
                    if f.contains(k) {
                        0.
                    } else {
                        weights[*k as usize] + SMALL_WEIGHT
                    }
                })
                .unwrap();
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
) -> Vec<Option<u8>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    if cluster_num <= 1 || data.len() <= 2 {
        return vec![None; data.len()];
    }
    assert_eq!(f.len(), data.len());
    let lim = limit / 3;
    let times = vec![0.005, 0.01, 0.10];
    let sum = data
        .iter()
        .flat_map(|e| e.id().bytes())
        .map(|e| e as u64)
        .sum::<u64>()
        + 1;
    let mut assignments = vec![];
    for (seed, pick_prob) in times.into_iter().enumerate() {
        let params = (lim, pick_prob, seed as u64 * sum);
        let (res, end) = gibbs_sampling_inner(data, labels, f, cluster_num, params, config, aln);
        assignments = res;
        if end {
            break;
        }
    }
    assignments
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
    let total_unit = data.iter().map(|e| e.len()).sum::<usize>();
    let (variants, pos) = select_variants(variants, total_unit, chain_len);
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
        let total_unit = data.iter().map(|e| e.len()).sum::<usize>();
        let (variants, pos) = select_variants(variants, total_unit, chain_len);
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
        if changed_num <= (data.len() as u32 / 30).max(3) {
            count += 1;
        } else {
            count = 0;
        }
        predictions.push_back(asn.clone());
        if predictions.len() as u32 > STABLE_LIMIT {
            predictions.pop_front();
        }
        report_gibbs(asn, answer, label, count, id, cluster_num, pick_prob);
        if (std::time::Instant::now() - start).as_secs() > limit {
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
    info!("Summary\t{}\t{}\t{:.3}\t{}", id, lp, beta, line.join("\t"));
}
pub fn soft_clustering_poa<F>(
    data: &[ERead],
    (label, answer): (&[u8], &[u8]),
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    config: &poa_hmm::Config,
    aln: &AlnParam<F>,
) -> Vec<Vec<f64>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    debug!("{}", cluster_num);
    if cluster_num <= 1 {
        return vec![vec![1.]; data.len()];
    }
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 99);
    let mut weights_of_reads =
        construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data);
    debug!("Model factory is built");
    let datasize = data.len() as f64;
    let border = label.len();
    let mut beta = INIT_BETA;
    for loop_num in 1.. {
        let init: Vec<Vec<_>> = (0..cluster_num)
            .map(|k| (0..k).map(|_| vec![0.; chain_len]).collect())
            .collect();
        let ws: Vec<f64> = (0..cluster_num)
            .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
            .collect();
        let wor = &weights_of_reads;
        let (variants, prev_lk): (Vec<_>, f64) = (0..REP_NUM)
            .map(|_| {
                let models: Vec<Vec<POA>> =
                    construct_models(wor, &data, cluster_num, chain_len, aln, &mut rng);
                variant_calling::variant_calling_all_pairs(&models, &data, config, &ws)
            })
            .fold((init, 0.), |(mut x, lk), (y, l)| {
                let div = REP_NUM as f64;
                x.iter_mut().zip(y).for_each(|(xss, yss)| {
                    xss.iter_mut().zip(yss).for_each(|(xs, ys)| {
                        xs.iter_mut().zip(ys).for_each(|(x, y)| *x += y * y / div);
                    })
                });
                (x, lk + l / REP_NUM as f64)
            });
        let betas = normalize_weights(&variants, beta);
        let wor = &mut weights_of_reads;
        update_weights(wor, border, &data, &betas, chain_len, config, &mut rng, aln);
        let lk = {
            let wor = &weights_of_reads;
            mf.update_seeds(&mut rng);
            let models: Vec<Vec<POA>> = (0..cluster_num)
                .map(|cl| mf.generate_model(&wor, &data, cl, aln))
                .collect();
            let ws: Vec<f64> = (0..cluster_num)
                .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
                .collect();
            variant_calling::get_lk(&models, &data, config, &ws)
        };
        info!("LK\t{}\t{}\t{}\t{:.3}", id, loop_num, lk, beta);
        report(id, &weights_of_reads, border, answer, cluster_num);
        let soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
        beta *= match lk.partial_cmp(&prev_lk).unwrap() {
            std::cmp::Ordering::Greater => BETA_INCREASE,
            _ if soe / datasize / (cluster_num as f64).ln() < ENTROPY_THR => break,
            _ => BETA_DECREASE,
        };
    }
    weights_of_reads
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
    _total: usize,
    chain: usize,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>) {
    let mut position = vec![false; chain];
    for bss in variants.iter_mut() {
        for bs in bss.iter_mut() {
            let thr = {
                let mut var = bs.clone();
                var.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let pos = (var.len() as f64 * 0.92).floor() as usize;
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
    // for (idx, bss) in variants.iter().enumerate() {
    //     for (jdx, bs) in bss.iter().enumerate() {
    //         let line: Vec<_> = bs
    //             .iter()
    //             .enumerate()
    //             .filter(|&x| *x.1 > 0.001)
    //             .map(|x| x.0)
    //             .map(|e| format!("{}", e))
    //             .collect();
    //         debug!("POS\t{}\t{}\t{}", idx, jdx, line.join(","));
    //     }
    // }
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

fn construct_models<F, R: Rng>(
    weights_of_reads: &[Vec<f64>],
    data: &[Read],
    cluster_num: usize,
    chain_len: usize,
    aln: &AlnParam<F>,
    rng: &mut R,
) -> Vec<Vec<POA>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let param = (aln.ins, aln.del, &aln.score);
    let mut chunks = vec![vec![vec![]; chain_len]; cluster_num];
    let clusters: Vec<_> = (0..cluster_num).collect();
    for _ in 0..data.len() {
        let i = rng.gen_range(0, data.len());
        let (read, weights) = (&data[i], &weights_of_reads[i]);
        let slots = *clusters.choose_weighted(rng, |&cl| weights[cl]).unwrap();
        for &(pos, ref unit) in read.iter() {
            chunks[slots][pos].push(unit.as_slice());
        }
    }
    let mut prior = vec![vec![]; chain_len];
    for read in data {
        if rng.gen_bool(PICK_PRIOR) {
            for &(pos, ref unit) in read.iter() {
                prior[pos].push(unit.as_slice());
            }
        }
    }
    let seeds: Vec<_> = rng.sample_iter(Standard).take(chain_len).collect();
    let prior: Vec<_> = prior
        .par_iter()
        .zip(seeds.par_iter())
        .map(|(cs, &s)| {
            let weights = vec![POA_PRIOR / cs.len() as f64; cs.len()];
            POA::default().update(cs, &weights, param, s)
        })
        .collect();
    let seeds: Vec<_> = rng.sample_iter(Standard).take(chain_len).collect();
    chunks
        .into_iter()
        .map(|css| {
            css.par_iter()
                .zip(prior.clone().into_par_iter())
                .zip(seeds.par_iter())
                .map(|((cs, poa), &seed)| poa.update(cs, &vec![1.; cs.len()], param, seed))
                .collect()
        })
        .collect()
}

fn likelihoods(models: &[Vec<POA>], read: &Read, c: &Config) -> Vec<Vec<(f64, usize)>> {
    models
        .iter()
        .map(|ms| {
            read.iter()
                .map(|&(pos, ref unit)| (ms[pos].forward(unit, c), pos))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn update_weights<F, R>(
    weights_of_reads: &mut [Vec<f64>],
    border: usize,
    data: &[Read],
    betas: &[Vec<Vec<f64>>],
    chain_len: usize,
    c: &Config,
    rng: &mut R,
    aln: &AlnParam<F>,
) where
    R: Rng,
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let cluster_num = weights_of_reads[0].len();
    let datasize = data.len() - border;
    let raw_likelihoods: Vec<_> =
        (0..REP_NUM).fold(vec![vec![vec![]; cluster_num]; datasize], |mut acc, _| {
            let models: Vec<Vec<POA>> =
                construct_models(weights_of_reads, &data, cluster_num, chain_len, aln, rng);
            data.par_iter()
                .skip(border)
                .zip(acc.par_iter_mut())
                .for_each(|(read, acc)| {
                    let likelihoods = likelihoods(&models, read, c);
                    for (cluster, lks) in likelihoods.into_iter().enumerate() {
                        acc[cluster].push(lks);
                    }
                });
            acc
        });
    let ws: Vec<f64> = (0..cluster_num)
        .map(|cl| weights_of_reads.iter().map(|e| e[cl]).sum::<f64>())
        .collect();
    weights_of_reads
        .par_iter_mut()
        .skip(border)
        .zip(raw_likelihoods.into_par_iter())
        .for_each(|(weights, likelihoods)| {
            weights.iter_mut().enumerate().for_each(|(k, w)| {
                *w = (0..cluster_num)
                    .map(|l| {
                        if k != l {
                            let (i, j) = (l.max(k), l.min(k));
                            let prior = digamma(ws[l]) - digamma(ws[k]);
                            let l_lk: Vec<_> = likelihoods[l]
                                .iter()
                                .map(|lks| {
                                    lks.iter().map(|&(l, p)| betas[i][j][p] * l).sum::<f64>()
                                })
                                .collect();
                            let l_lk = logsumexp(&l_lk);
                            let k_lk: Vec<_> = likelihoods[k]
                                .iter()
                                .map(|lks| {
                                    lks.iter().map(|&(l, p)| betas[i][j][p] * l).sum::<f64>()
                                })
                                .collect();
                            let k_lk = logsumexp(&k_lk);
                            prior + l_lk - k_lk
                        } else {
                            0.
                        }
                    })
                    .map(|lkdiff| lkdiff.exp())
                    .sum::<f64>()
                    .recip();
            });
            let sum = weights.iter().sum::<f64>();
            weights.iter_mut().for_each(|e| *e /= sum);
        });
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

/// Return likelihood of the assignments.
pub fn bic<F>(
    data: &[ERead],
    weights_of_reads: &[Vec<f64>],
    cluster_num: usize,
    config: &Config,
    aln: &AlnParam<F>,
) -> (f64, usize, usize)
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
        .collect();
    let lk = likelihood_of_models(&models, &data, &ws, config);
    (lk, (cluster_num - 1) * (1 + data.len()), data.len())
}

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
    bic(data, weights_of_reads, cluster_num, config, aln).0
}

fn likelihood_of_models(models: &[Vec<POA>], data: &[Read], ws: &[f64], c: &Config) -> f64 {
    let len = ws.iter().sum::<f64>();
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
                    lk + digamma(*w) - digamma(len)
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

pub fn diff_between<F>(
    data: &[ERead],
    wor: &[Vec<f64>],
    cluster_num: usize,
    aln: &AlnParam<F>,
) -> Vec<Vec<i32>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 9999);
    let models: Vec<Vec<POA>> = construct_models(wor, &data, cluster_num, chain_len, aln, &mut rng);
    (0..cluster_num)
        .map(|i| {
            (0..cluster_num)
                .map(|j| {
                    models[i]
                        .iter()
                        .zip(&models[j])
                        .map(|(m_i, m_j)| local_alignment(&m_i.consensus(), &m_j.consensus()))
                        .sum::<i32>()
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn local_alignment(xs: &[u8], ys: &[u8]) -> i32 {
    let (ins, del, score) = (-1, -1, |x, y| if x == y { 1 } else { -1 });
    let mut prev = vec![0; xs.len() + 1];
    let mut updated = vec![0; xs.len() + 1];
    let mut max_so_far = 0;
    for &y in ys {
        for (j, &x) in xs.iter().enumerate() {
            updated[j + 1] = (prev[j] + score(x, y))
                .max(prev[j + 1] + ins)
                .max(updated[j] + del)
                .max(0);
        }
        max_so_far = updated.iter().copied().max().unwrap().max(max_so_far);
        std::mem::swap(&mut prev, &mut updated);
    }
    max_so_far
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
        let mut choises = vec![true; cluster_num];
        let forbidden: &Vec<u8> = &forbidden[idx + border];
        forbidden
            .iter()
            .for_each(|&cl| choises[cl as usize] = false);
        let choises: Vec<_> = choises
            .into_iter()
            .enumerate()
            .filter_map(|(idx, b)| if b { Some(idx) } else { None })
            .collect();
        let mut bucket = vec![0; cluster_num];
        (0..num_of_ball).for_each(|_| bucket[*choises.choose(&mut rng).unwrap()] += 1);
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
    let total_unit = data.iter().map(|e| e.len()).sum::<usize>();
    let (variants, pos) = select_variants(variants, total_unit, chain_len);
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
