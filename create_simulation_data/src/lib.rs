extern crate dbg_hmm;
extern crate last_tiling;
extern crate log;
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

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    let (x1, x2) = ((x1 - log_denominator).exp(), (x2 - log_denominator).exp());
    assert!((x1 + x2 - 1.).abs() < 0.0001, "{}", x1 + x2);
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
