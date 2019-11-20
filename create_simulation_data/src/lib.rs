extern crate dbg_hmm;
extern crate last_tiling;
extern crate log;
extern crate rayon;
use dbg_hmm::*;
use last_tiling::EncodedRead;

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

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    let (x1, x2) = ((x1 - log_denominator).exp(), (x2 - log_denominator).exp());
    assert!((x1 + x2 - 1.).abs() < 0.0001, "{}", x1 + x2);
    (x1, x2)
}

pub type Seq = Vec<u8>;
// UnitID -> Sequences -> (color, seq)
pub fn construct_predictor(
    training: &[EncodedRead],
    is_original: &[bool],
    len: usize,
) -> Vec<(Vec<Seq>, Vec<Seq>)> {
    let mut chunks = vec![(vec![], vec![]); len];
    use last_tiling::revcmp;
    for (&is_original, read) in is_original.iter().zip(training.iter()) {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                revcmp(unit.bases.as_bytes())
            };
            if is_original {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    chunks
}

pub fn unit_predict(query: &[u8], templates: &[Seq], k: usize) -> f64 {
    DBGHMM::new(templates, k).forward(&query, &DEFAULT_CONFIG)
}

pub fn unit_predict_by(query: &[u8], refs: &[Seq], k: usize, c: &Config) -> f64 {
    DBGHMM::new(refs, k).forward(&query, c)
}
