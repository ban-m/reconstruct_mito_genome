extern crate dbg_hmm;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate rand_xoshiro;
use rand_xoshiro::Xoshiro256StarStar;
extern crate rand;
use rand::{Rng, SeedableRng};
extern crate rayon;
use dbg_hmm::*;
use last_tiling::EncodedRead;

pub const BADREAD_CONFIG: dbg_hmm::Config = dbg_hmm::Config {
    mismatch: 0.0344,
    p_match: 0.88,
    p_ins: 0.0549,
    p_del: 0.0651,
    p_extend_ins: 0.0337,
    p_extend_del: 0.1787,
    p_del_to_ins: 0.0,
    match_score: 1,
    mism_score: -1,
    del_score: -1,
    ins_score: -1,
    base_freq: [0.25, 0.25, 0.25, 0.25],
};

/// A simple repr for EncodedRead.
#[derive(Debug, Clone)]
pub struct ERead {
    pub id: String,
    pub seq: Vec<CUnit>,
}

#[derive(Debug, Clone)]
pub struct CUnit {
    pub contig: u16,
    pub unit: u16,
    pub bases: Vec<u8>,
}

impl CUnit {
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
}

impl ERead {
    pub fn new(er: EncodedRead) -> Self {
        let EncodedRead { id, seq } = er;
        let seq = seq
            .iter()
            .filter_map(|u| u.encode())
            .map(|e| {
                let contig = e.contig;
                let unit = e.unit;
                let bases = if e.is_forward {
                    e.bases.as_bytes().to_vec()
                } else {
                    last_tiling::revcmp(e.bases.as_bytes())
                };
                CUnit {
                    contig,
                    unit,
                    bases,
                }
            })
            .collect();
        Self { id, seq }
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn seq(&self) -> &[CUnit] {
        &self.seq
    }
}

pub fn predict_by_naive(ps: &[(f64, f64)]) -> u8 {
    let p1_larget_than_p2 = ps
        .iter()
        .map(|(l1, l2)| l1 - l2)
        .sum::<f64>()
        .is_sign_positive();
    if p1_larget_than_p2 {
        0
    } else {
        1
    }
}

pub fn predict_by_independent(ps: &[(f64, f64)]) -> u8 {
    let sum = ps.iter().filter(|(l1, l2)| l1 > l2).count();
    if sum * 2 >= ps.len() {
        0
    } else {
        1
    }
}

pub fn predict_by_sol(ps: &[(f64, f64)]) -> u8 {
    let xlnx = |x: &f64| if x == &0. { 0. } else { x * x.ln() };
    let entropy: Vec<_> = ps
        .iter()
        .map(|(o, m)| {
            let (o, m) = as_weight(*o, *m);
            (2f64).ln() + xlnx(&o) + xlnx(&m)
        })
        .collect();
    let p1_larger_than_p2 = ps
        .iter()
        .zip(entropy.iter())
        .map(|((l1, l2), w)| w * (l1 - l2))
        .sum::<f64>()
        .is_sign_positive();
    if p1_larger_than_p2 {
        0
    } else {
        1
    }
}

pub fn predict_by_sow(ps: &[(f64, f64)]) -> u8 {
    let xlnx = |x: &f64| if x == &0. { 0. } else { x * x.ln() };
    let entropy: Vec<_> = ps
        .iter()
        .map(|(o, m)| {
            let (o, m) = as_weight(*o, *m);
            (2f64).ln() + xlnx(&o) + xlnx(&m)
        })
        .collect();
    let p1_lager_than_p2 = ps
        .iter()
        .map(|(l1, l2)| as_weight(*l1, *l2))
        .zip(entropy.iter())
        .map(|((w1, w2), w)| (w1 - w2) * w)
        .sum::<f64>()
        .is_sign_positive();
    if p1_lager_than_p2 {
        0
    } else {
        1
    }
}

pub fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    let (x1, x2) = ((x1 - log_denominator).exp(), (x2 - log_denominator).exp());
    (x1, x2)
}

pub type Seq = Vec<u8>;
// UnitID -> Sequences -> (color, seq)
pub fn construct_predictor<'a>(
    training: &'a [ERead],
    is_original: &[bool],
    len: usize,
) -> Vec<(Vec<&'a [u8]>, Vec<&'a [u8]>)> {
    let mut chunks = vec![(vec![], vec![]); len];
    for (&is_original, read) in is_original.iter().zip(training.iter()) {
        for unit in read.seq.iter() {
            let u = unit.unit as usize;
            let seq = unit.bases.as_slice();
            if is_original {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    chunks
}

pub fn unit_predict(query: &[u8], templates: &[&[u8]], k: usize) -> f64 {
    DBGHMM::new_from_ref(templates, k).forward(&query, &DEFAULT_CONFIG)
}

pub fn unit_predict_by(query: &[u8], refs: &[&[u8]], k: usize, c: &Config) -> f64 {
    DBGHMM::new_from_ref(refs, k).forward(&query, c)
}

/// Conpute Log (w0 * exp(l0) + w1 * exp(l1)). If you want to do the same thing for
/// two vectors, see calc_logsum_vec
pub fn calc_logsum(l0: f64, l1: f64, w0: f64, w1: f64) -> f64 {
    // Log (w0 * exp(l0) + w1 * exp(l1))
    let f1 = w0.ln() + l0;
    let f2 = w1.ln() + l1;
    let m = f1.max(f2);
    m + ((f1 - m).exp() + (f2 - m).exp()).ln()
}

/// Log Sum w_i * exp(l_i). For two pairs of values, use calc_logsum instead.
pub fn calc_logsum_vec(ls: &[f64], ws: &[f64]) -> f64 {
    let m = ls
        .iter()
        .zip(ws.iter())
        .map(|(&l, &w)| l + w.ln())
        .fold(std::f64::MIN, |x, y| x.min(y));
    m + ls
        .iter()
        .zip(ws.iter())
        .map(|(&l, &w)| l + w.ln() - m)
        .map(|e| e.exp())
        .sum::<f64>()
        .ln()
}

/// Generate dataset. Return Reads(Chunked), assignment(label for training data),
/// answer(label for test data), and l(the length of training data).
/// The data is ordered from training set followed by test set.
pub fn generate_dataset<T: Rng>(
    template1: &[Vec<u8>],
    template2: &[Vec<u8>],
    coverage: usize,
    test_num: usize,
    rng: &mut T,
    skew: f64,
) -> (Vec<Vec<Vec<u8>>>, Vec<u8>, Vec<u8>, usize) {
    let answer: Vec<_> = (0..test_num + coverage)
        .map(|i| {
            if i < 4 {
                i % 2 == 0
            } else {
                rng.gen_bool(skew)
            }
        })
        .map(|b| if b { 0 } else { 1 })
        .collect();
    let mut gen = |t: &[Vec<u8>]| {
        t.iter()
            .map(|e| gen_sample::introduce_randomness(e, rng, &gen_sample::PROFILE))
            .collect::<Vec<_>>()
    };
    debug!("Coverage:{}\ttest_num:{}", coverage, test_num);
    let dataset: Vec<_> = answer
        .iter()
        .map(|e| {
            if e % 2 == 0 {
                gen(template1)
            } else {
                gen(template2)
            }
        })
        .collect();
    let assignment: Vec<_> = answer.iter().take(coverage).copied().collect();
    let answer: Vec<_> = answer.into_iter().skip(coverage).collect();
    assert_eq!(dataset.len(), assignment.len() + answer.len());
    let l = assignment.len();
    (dataset, assignment, answer, l)
}

pub fn naive_solve(dataset: &[Vec<Vec<u8>>], label: &[u8], border: usize, k: usize) -> Vec<u8> {
    let (d1, d2): (Vec<_>, Vec<_>) = dataset
        .iter()
        .take(border)
        .zip(label)
        .partition(|&(_, &ans)| ans == 0);
    let mut d1: Vec<_> = d1.into_iter().map(|(x, _)| x).collect();
    let mut d2: Vec<_> = d2.into_iter().map(|(x, _)| x).collect();
    let mut model1 = construct_from_reads(&d1, k);
    let mut model2 = construct_from_reads(&d2, k);
    let mut pred = vec![0; dataset.len() - border];
    let mut is_updated = true;
    while is_updated {
        is_updated = pred
            .iter_mut()
            .enumerate()
            .map(|(idx, p)| {
                let lkdiff = dataset[idx + border]
                    .iter()
                    .enumerate()
                    .map(|(idx, chunk)| {
                        model1[idx].forward(chunk, &DEFAULT_CONFIG)
                            - model2[idx].forward(chunk, &DEFAULT_CONFIG)
                    })
                    .sum::<f64>()
                    .is_sign_positive();
                let up = if lkdiff { 0 } else { 1 };
                let res = *p != up;
                *p = up;
                res
            })
            .fold(false, |p, q| p | q);
        d1 = dataset
            .iter()
            .zip(label.iter().chain(pred.iter()))
            .filter(|&(_, &b)| b == 0)
            .map(|(e, _)| e)
            .collect();
        d2 = dataset
            .iter()
            .zip(label.iter().chain(pred.iter()))
            .filter(|&(_, &b)| b == 1)
            .map(|(e, _)| e)
            .collect();
        model1 = construct_from_reads(&d1, k);
        model2 = construct_from_reads(&d2, k);
    }
    pred
}

/// Predict by EM algorithm. Return the length of return value is the number of test case.
pub fn em_solve(
    dataset: &[Vec<Vec<u8>>],
    label: &[u8],
    border: usize,
    k: usize,
    s: u64,
) -> Vec<u8> {
    let pred = naive_solve(dataset, label, border, k);
    let (mut w0, mut w1) = {
        let (w0, w1) = label.iter().chain(pred.iter()).fold((0, 0), |(x, y), &p| {
            if p == 0 {
                (x + 1, y)
            } else {
                (x, y + 1)
            }
        });
        let len = dataset.len() as f64;
        (w0 as f64 / len, w1 as f64 / len)
    };
    let mut gamma0: Vec<_> = label
        .iter()
        .chain(pred.iter())
        .map(|&p| if p == 0 { 1. } else { 0. })
        .collect();
    let mut gamma1: Vec<_> = label
        .iter()
        .chain(pred.iter())
        .map(|&p| if p == 0 { 0. } else { 1. })
        .collect();
    let mut model0 = construct_with_weights(dataset, &gamma0, k);
    let mut model1 = construct_with_weights(dataset, &gamma1, k);
    let mut lk = calc_lk(&dataset[border..], &model0, &model1, w0, w1);
    info!("LK\t{:.3}\t0", lk);
    //    use std::io::{BufWriter, Write};
    let ep = 0.0001;
    for i in 0..10 {
        for (idx, read) in dataset.iter().enumerate().skip(border) {
            let log_m0 = read
                .iter()
                .enumerate()
                .map(|(idx, chunk)| model0[idx].forward(chunk, &DEFAULT_CONFIG))
                .sum::<f64>();
            let log_m1 = read
                .iter()
                .enumerate()
                .map(|(idx, chunk)| model1[idx].forward(chunk, &DEFAULT_CONFIG))
                .sum::<f64>();
            let w = calc_logsum(log_m0, log_m1, w0, w1);
            gamma0[idx] = (w0.ln() + log_m0 - w).exp();
            gamma1[idx] = (w1.ln() + log_m1 - w).exp();
        }
        let tot = gamma0.iter().sum::<f64>() + gamma1.iter().sum::<f64>();
        w0 = gamma0.iter().sum::<f64>() / tot;
        w1 = gamma1.iter().sum::<f64>() / tot;
        model0 = construct_with_weights(&dataset, &gamma0, k);
        model1 = construct_with_weights(&dataset, &gamma1, k);
        let next_lk = calc_lk(&dataset[border..], &model0, &model1, w0, w1);
        info!("LK\t{:.3}\t{}", next_lk, i);
        if (next_lk - lk).abs() < ep {
            break;
        } else {
            lk = next_lk;
        }
    }
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(s + 9893242349);
    gamma0
        .iter()
        .skip(border)
        .map(|&g0| if rng.gen_range(0., 1.) < g0 { 0 } else { 1 })
        .collect()
}

pub fn calc_lk(ds: &[Vec<Vec<u8>>], m0: &[DBGHMM], m1: &[DBGHMM], w0: f64, w1: f64) -> f64 {
    ds.iter()
        .map(|read| {
            read.iter()
                .enumerate()
                .map(|(idx, chunk)| {
                    let m0 = m0[idx].forward(chunk, &DEFAULT_CONFIG);
                    let m1 = m1[idx].forward(chunk, &DEFAULT_CONFIG);
                    calc_logsum(m0, m1, w0, w1)
                })
                .sum::<f64>()
        })
        .sum::<f64>()
}

pub fn construct_from_reads(ds: &[&Vec<Vec<u8>>], k: usize) -> Vec<DBGHMM> {
    let len = ds[0].len();
    let mut res: Vec<Vec<&[u8]>> = vec![vec![]; len];
    assert!(ds.iter().all(|e| e.len() == len));
    for read in ds {
        for (idx, chunk) in read.iter().enumerate() {
            res[idx].push(chunk);
        }
    }
    let mut f = Factory::new();
    res.into_iter()
        .map(|e| f.generate_from_ref(&e, k))
        .collect()
}

pub fn construct_with_weights(ds: &[Vec<Vec<u8>>], ws: &[f64], k: usize) -> Vec<DBGHMM> {
    let len = ds[0].len();

    let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; len];
    for read in ds.into_iter() {
        for (idx, chunk) in read.iter().enumerate() {
            chunks[idx].push(chunk);
        }
    }
    let mut f = Factory::new();
    chunks
        .into_iter()
        .map(|cs| f.generate_with_weight(&cs, &ws, k))
        .collect()
}
