//! This module is to call variant model(i.e. the positions where some models have variant site).
//! The method is based on diagonalization of total distance between the center axis and
//! each datapoint.
use super::ERead;
use dbg_hmm::Config;
use dbg_hmm::DBGHMM;
use na::DMatrix;
use rayon::prelude::*;
use std::time::Instant;
/// return the weight of each position.
pub fn variant_call(
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    c: &Config,
    pos: &[Vec<usize>],
    num_chain: usize,
) -> Vec<Vec<f64>> {
    let s1 = Instant::now();
    let matrix = calc_matrix(models, data, c, &pos, num_chain);
    let s2 = Instant::now();
    let eigens = matrix.clone().symmetric_eigen();
    let s3 = Instant::now();
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    eprintln!("{:?}", max);
    debug!("LK:Eigen={:?}:{:?}", s2 - s1, s3 - s2);
    let betas: Vec<_> = eigens
        .eigenvectors
        .column(max.0)
        .iter()
        .map(|e| e * e)
        .collect();
    pos.iter()
        .map(|e| e.iter().map(|&idx| betas[idx]).collect())
        .collect()
}

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
                res[contig][i] = len + i;
            }
            len += max - min + 1;
        }
    }
    (res, len)
}

fn calc_matrix(
    models: &[Vec<Vec<DBGHMM>>],
    data: &[ERead],
    c: &Config,
    pos: &[Vec<usize>],
    chain_len: usize,
) -> DMatrix<f64> {
    let num_cluster = models.len();
    data.par_iter()
        .map(|read| lks(models, read, c, pos, chain_len))
        .map(|data| DMatrix::from_row_slice(num_cluster, chain_len, &data))
        .fold(
            || DMatrix::zeros(chain_len, chain_len),
            |x, l| {
                let trans = l.clone().transpose();
                let ltl = trans.clone() * l.clone();
                let ones = DMatrix::repeat(num_cluster, num_cluster, 1.);
                let lt_one_l = trans * ones * l;
                x + ltl - lt_one_l / num_cluster as f64
            },
        )
        .reduce(|| DMatrix::zeros(chain_len, chain_len), |x, y| x + y)
}

fn lks(
    models: &[Vec<Vec<DBGHMM>>],
    read: &ERead,
    c: &Config,
    pos: &[Vec<usize>],
    chain_len: usize,
) -> Vec<f64> {
    models
        .iter()
        .flat_map(|m| lk(m, read, c, pos, chain_len))
        .collect()
}

fn lk(
    m: &[Vec<DBGHMM>],
    read: &ERead,
    c: &Config,
    pos: &[Vec<usize>],
    chain_len: usize,
) -> Vec<f64> {
    let mut res: Vec<_> = vec![0.; chain_len];
    for unit in read.seq.iter() {
        let pos = pos[unit.contig()][unit.unit()];
        let lk = m[unit.contig()][unit.unit()].forward(unit.bases(), c);
        assert!(res[pos] == 0.);
        res[pos] = lk;
    }
    res
}
