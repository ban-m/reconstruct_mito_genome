extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate poa_hmm;
extern crate rand_xoshiro;
use poa_hmm::gen_sample;
extern crate rand;
use rand::Rng;
extern crate rayon;
use last_tiling::EncodedRead;
use rayon::prelude::*;

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
        let EncodedRead { id, seq, .. } = er;
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
    let xlnx = |x: &f64| {
        if x.abs() < f64::EPSILON {
            0.
        } else {
            x * x.ln()
        }
    };
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
    let xlnx = |x: &f64| {
        if x.abs() < f64::EPSILON {
            0.
        } else {
            x * x.ln()
        }
    };
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

// pub fn construct_predictor<'a>(
//     training: &'a [ERead],
//     is_original: &[bool],
//     len: usize,
// ) -> Vec<(Vec<&'a [u8]>, Vec<&'a [u8]>)> {
//     let mut chunks = vec![(vec![], vec![]); len];
//     for (&is_original, read) in is_original.iter().zip(training.iter()) {
//         for unit in read.seq.iter() {
//             let u = unit.unit as usize;
//             let seq = unit.bases.as_slice();
//             if is_original {
//                 chunks[u].0.push(seq);
//             } else {
//                 chunks[u].1.push(seq);
//             }
//         }
//     }
//     chunks
// }

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
pub fn generate_mul_data<T: Rng>(
    templates: &[Vec<Vec<u8>>],
    coverage: usize,
    test_num: usize,
    rng: &mut T,
    probs: &[f64],
    profile: &gen_sample::Profile,
) -> (Vec<Vec<Vec<u8>>>, Vec<u8>, Vec<u8>, usize) {
    let answer: Vec<_> = (0..probs.len()).flat_map(|i| vec![i; coverage]).collect();
    let answer: Vec<_> = answer
        .into_iter()
        .chain(probs.iter().enumerate().flat_map(|(idx, &prob)| {
            let num = (test_num as f64 * prob).ceil() as usize;
            vec![idx; num]
        }))
        .map(|e| e as u8)
        .collect();
    let mut gen = |t: &[Vec<u8>]| {
        t.iter()
            .map(|e| gen_sample::introduce_randomness(e, rng, profile))
            .collect::<Vec<_>>()
    };
    debug!("Coverage:{}\tTest Num:{}", coverage, test_num);
    let dataset: Vec<_> = answer
        .iter()
        .map(|&idx| gen(&templates[idx as usize]))
        .collect();
    let border = probs.len() * coverage;
    let (label, answer) = answer.split_at(border);
    assert_eq!(dataset.len(), label.len() + answer.len());
    (dataset, label.to_vec(), answer.to_vec(), border)
}

/// Generate dataset. Return Reads(Chunked), assignment(label for training data),
/// answer(label for test data), and l(the length of training data).
/// The data is ordered from training set followed by test set.
pub fn generate_dataset<T: Rng>(
    template0: &[Vec<u8>],
    template1: &[Vec<u8>],
    coverage: usize,
    test_num: usize,
    rng: &mut T,
    skew: f64,
) -> (Vec<Vec<Vec<u8>>>, Vec<u8>, Vec<u8>, usize) {
    let answer: Vec<_> = (0..test_num + coverage)
        .map(|i| {
            if i < 6 {
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
        .map(|&e| {
            if e == 0 {
                gen(template0)
            } else {
                gen(template1)
            }
        })
        .collect();
    let assignment: Vec<_> = answer.iter().take(coverage).copied().collect();
    let answer: Vec<_> = answer.into_iter().skip(coverage).collect();
    assert_eq!(dataset.len(), assignment.len() + answer.len());
    let l = assignment.len();
    (dataset, assignment, answer, l)
}

pub fn align_solve(data: &[Vec<Vec<u8>>], label: &[u8], border: usize) -> Vec<u8> {
    let data = data
        .iter()
        .map(|read| {
            read.iter()
                .flat_map(|e| e.iter().map(|&e| e))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut d0: Vec<&Vec<_>> = data
        .iter()
        .zip(label.iter())
        .filter_map(|(r, &b)| if b == 0 { Some(r) } else { None })
        .collect();
    let mut d1: Vec<&Vec<_>> = data
        .iter()
        .zip(label.iter())
        .filter_map(|(r, &b)| if b != 0 { Some(r) } else { None })
        .collect();
    let mut assign: Vec<_> = vec![1; data.len()];
    let mut updated = true;
    while updated {
        use bio_utils::alignments::edit_dist;
        updated = assign
            .par_iter_mut()
            .enumerate()
            .skip(border)
            .map(|(idx, a)| {
                let query = &data[idx];
                let min0 = d0.iter().map(|r| edit_dist(r, query)).min();
                let min1 = d1.iter().map(|r| edit_dist(r, query)).min();
                let p = if min0 < min1 { 0 } else { 1 };
                let updated = p != *a;
                *a = p;
                updated
            })
            .reduce(|| false, |p, q| q | p);
        //.fold(false, |p, q| q | p);
        d0 = data
            .iter()
            .zip(assign.iter())
            .filter_map(|(r, &b)| if b == 0 { Some(r) } else { None })
            .collect();
        d1 = data
            .iter()
            .zip(assign.iter())
            .filter_map(|(r, &b)| if b != 0 { Some(r) } else { None })
            .collect();
    }
    assign.into_iter().skip(border).collect()
}

/// Return log ( x.exp() + y.exp())
pub fn logsumexp(x: f64, y: f64) -> f64 {
    let max = x.max(y);
    max + ((x - max).exp() + (y - max).exp()).ln()
}
