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
    centrize: bool,
) -> Vec<f64> {
    let matrix = calc_matrix_poa(models, data, c, centrize);
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

// Return likelihoods of the read.
// the cl * i + j -th element is the likelihood for the i-th cluster at the j-th position.
fn lks_poa(models: &[Vec<POA>], read: &super::Read, c: &poa_hmm::Config, cl: usize) -> Vec<f64> {
    let mut res = vec![0.; cl * models.len()];
    for (i, ms) in models.iter().enumerate() {
        for &(pos, ref unit) in read.iter() {
            res[i * cl + pos] = ms[pos].forward(unit, c);
        }
    }
    res
}

// Return the centrize vector for matrices.
// In other words, it calculates arrays of vector p_1,...,p_{row}, where
// each c_i is column-dimensional vector.
fn centrize_vector_of(
    data: &[super::Read],
    matrices: &[Vec<f64>],
    row: usize,
    column: usize,
) -> Vec<Vec<f64>> {
    let mut centrize_vector = vec![vec![0.; column]; row];
    // How many units is considered at position i.
    let mut counts = vec![0; column];
    // Compute the avarage vector at each **column(position)**
    for (read, matrix) in data.iter().zip(matrices.iter()) {
        for &(pos, _) in read.iter() {
            counts[pos] += 1;
            for cluster in 0..row {
                centrize_vector[cluster][pos] += matrix[column * cluster + pos];
            }
        }
    }
    centrize_vector.iter_mut().for_each(|row_vec| {
        row_vec
            .iter_mut()
            .zip(counts.iter())
            .for_each(|(sum, &count)| *sum = *sum / count as f64)
    });
    // Compute projection to axis, i.e., <c,1>/K 1 - c.
    // Here, c is the column vector at each position.
    // This is a centrize vector.
    // First, compute the <c,1> for each column.
    // We could have this by just add by columnwise.
    let mut inner_product = centrize_vector
        .iter()
        .fold(vec![0.; column], |mut acc, row_vec| {
            acc.iter_mut()
                .zip(row_vec.iter())
                .for_each(|(x, y)| *x += y);
            acc
        });
    // Then, normalize it by dividing by length of row(number of cluster).
    inner_product.iter_mut().for_each(|x| *x /= row as f64);
    // Lastly, compute <c,1>/K 1 - c for each position.
    centrize_vector.iter_mut().for_each(|cv| {
        cv.iter_mut()
            .zip(inner_product.iter())
            .for_each(|(x, prod)| *x = prod - *x);
    });
    centrize_vector
}

fn calc_matrix_poa(
    models: &[Vec<POA>],
    data: &[super::Read],
    c: &poa_hmm::Config,
    centrize: bool,
) -> DMatrix<f64> {
    let num_cluster = models.len();
    let chain_len = models[0].len();
    // LK matrix of each read, i.e., lk_matrix[i] would
    // be lk matrix for the i-th read.
    // Note that the matrix is flattened.
    // To access the likelihood of the j-th position of the k-th cluster,
    // lk_matrix[i][chan_len * k + j] would work.
    let mut lk_matrices: Vec<Vec<f64>> = data
        .par_iter()
        .map(|read| lks_poa(models, read, c, chain_len))
        .collect();
    if centrize {
        let cv = centrize_vector_of(data, &lk_matrices, num_cluster, chain_len);
        // Centrize the matrices. In other words,
        // we add the offset cv for each position.
        // Note that we only care at the position where the likelihood is non-zero value.
        lk_matrices
            .par_iter_mut()
            .zip(data.par_iter())
            .for_each(|(matrix, read)| {
                for &(pos, _) in read.iter() {
                    for row in 0..num_cluster {
                        matrix[row * chain_len + pos] += cv[row][pos]
                    }
                }
            });
    };
    lk_matrices
        .into_par_iter()
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
