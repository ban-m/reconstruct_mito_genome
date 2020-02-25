//! This module is to call variant model(i.e. the positions where some models have variant site).
//! The method is based on diagonalization of total distance between the center axis and
//! each datapoint.
use dbg_hmm::Config;
use dbg_hmm::DBGHMM;
use na::DMatrix;
use poa_hmm::POA;
use rayon::prelude::*;
/// return the weight of each position.
pub fn variant_call(models: &[Vec<DBGHMM>], data: &[super::Read], c: &Config) -> Vec<f64> {
    let matrix = calc_matrix(models, data, c);
    let eigens = matrix.symmetric_eigen();
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    eigens.eigenvectors.column(max.0).iter().copied().collect()
}

/// return the weight of each position.
pub fn variant_call_poa(
    models: &[Vec<POA>],
    data: &[super::Read],
    c: &poa_hmm::Config,
) -> Vec<f64> {
    let matrix = calc_matrix_poa(models, data, c);
    let eigens = matrix.symmetric_eigen();
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    eigens.eigenvectors.column(max.0).iter().copied().collect()
}

fn calc_matrix(models: &[Vec<DBGHMM>], data: &[super::Read], c: &Config) -> DMatrix<f64> {
    let num_cluster = models.len();
    let chain_len = models[0].len();
    data.par_iter()
        .map(|read| lks(models, read, c, chain_len))
        .map(|data| DMatrix::from_row_slice(num_cluster, chain_len, &data))
        .fold(
            || DMatrix::zeros(chain_len, chain_len),
            |x, l| {
                let trans = l.clone().transpose();
                let ones = DMatrix::repeat(num_cluster, num_cluster, 1.);
                let unit = DMatrix::identity(num_cluster, num_cluster);
                let reg = unit - ones / num_cluster as f64;
                x + trans * reg * l
            },
        )
        .reduce(|| DMatrix::zeros(chain_len, chain_len), |x, y| x + y)
}

fn lks(models: &[Vec<DBGHMM>], read: &super::Read, c: &Config, cl: usize) -> Vec<f64> {
    let mut data: Vec<_> = models.iter().flat_map(|m| lk(m, read, c, cl)).collect();
    let ng_column: Vec<_> = data.chunks_exact(cl).flat_map(check_lk).collect();
    let class = models.len();
    for column in ng_column {
        for i in 0..class {
            data[i * cl + column] = 0.;
        }
    }
    data
}

fn check_lk(data: &[f64]) -> Vec<usize> {
    use super::LK_LIMIT;
    data.iter()
        .enumerate()
        .filter_map(|(idx, &lk)| if lk < LK_LIMIT { Some(idx) } else { None })
        .collect()
}

fn lk(m: &[DBGHMM], read: &super::Read, c: &Config, cl: usize) -> Vec<f64> {
    let mut res: Vec<_> = vec![0.; cl];
    for &(pos, ref unit) in read.iter() {
        res[pos] = m[pos].forward(unit, c);
    }
    res
}

fn calc_matrix_poa(models: &[Vec<POA>], data: &[super::Read], c: &poa_hmm::Config) -> DMatrix<f64> {
    let num_cluster = models.len();
    let chain_len = models[0].len();
    let lk_poa = |m: &[POA], read: &super::Read, c, cl| {
        let mut res: Vec<_> = vec![0.; cl];
        for &(pos, ref unit) in read.iter() {
            res[pos] = m[pos].forward(unit, c);
        }
        res
    };

    let lks_poa = |models: &[Vec<POA>], read, c, cl| {
        let mut data: Vec<_> = models.iter().flat_map(|m| lk_poa(m, read, c, cl)).collect();
        let ng_column: Vec<_> = data.chunks_exact(cl).flat_map(check_lk).collect();
        let class = models.len();
        for column in ng_column {
            for i in 0..class {
                data[i * cl + column] = 0.;
            }
        }
        data
    };
    data.par_iter()
        .map(|read| lks_poa(models, read, c, chain_len))
        .map(|data| DMatrix::from_row_slice(num_cluster, chain_len, &data))
        .fold(
            || DMatrix::zeros(chain_len, chain_len),
            |x, l| {
                let trans = l.clone().transpose();
                let ones = DMatrix::repeat(num_cluster, num_cluster, 1.);
                let unit = DMatrix::identity(num_cluster, num_cluster);
                let reg = unit - ones / num_cluster as f64;
                x + trans * reg * l
            },
        )
        .reduce(|| DMatrix::zeros(chain_len, chain_len), |x, y| x + y)
}
