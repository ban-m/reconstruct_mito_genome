//! This module is to call variant model(i.e. the positions where some models have variant site).
//! The method is based on diagonalization of total distance between the center axis and
//! each datapoint.
use na::DMatrix;
use poa_hmm::POA;
use rayon::prelude::*;

fn calc_matrix_poa(
    models: &[Vec<POA>],
    data: &[super::Read],
    to_centrize: bool,
    ws: &[f64],
    c: &super::poa_clustering::ClusteringConfig,
) -> (Vec<Vec<f64>>, f64) {
    let (lk_matrices, lk) = if to_centrize {
        let (matrices, lk) = calc_matrices_poa(models, data, ws, c);
        (centrize(matrices, c.cluster_num, c.chain_len), lk)
    } else {
        calc_matrices_poa(models, data, ws, c)
    };
    (lk_matrices, lk)
}

// Calculate LK matrices for each read. Total LK would be also returned.
// Note that each matrix is flattened in row-order.
fn calc_matrices_poa(
    models: &[Vec<POA>],
    data: &[super::Read],
    ws: &[f64],
    c: &super::poa_clustering::ClusteringConfig,
) -> (Vec<Vec<f64>>, f64) {
    // LK matrix of each read, i.e., lk_matrix[i] would
    // be lk matrix for the i-th read.
    // Note that the matrix is flattened.
    // To access the likelihood of the j-th position of the k-th cluster,
    // lk_matrix[i][chain_len * k + j] would work.
    // !!!!!!!!!!!
    let read_to_lks = |read| lks_poa(models, read, &c.poa_config, c.chain_len);
    let lk_matrices: Vec<Vec<f64>> = if c.is_par {
        data.par_iter().map(read_to_lks).collect()
    } else {
        data.iter().map(read_to_lks).collect()
    };
    let lk = lk_matrices
        .iter()
        .map(|matrix| {
            assert_eq!(matrix.len() / c.chain_len, c.cluster_num);
            let lks: Vec<_> = matrix
                .chunks_exact(c.chain_len)
                .zip(ws.iter())
                .map(|(chunks, w)| w.ln() + chunks.iter().sum::<f64>())
                .collect();
            crate::utils::logsumexp(&lks)
        })
        .sum::<f64>();
    (lk_matrices, lk)
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

fn centrize(mut matrices: Vec<Vec<f64>>, row: usize, column: usize) -> Vec<Vec<f64>> {
    let cv = centrize_vector_of(&matrices, row, column);
    // Centrize the matrices. In other words,
    // we add the offset cv for each position.
    // Note that we only care at the position where the likelihood is non-zero value.
    matrices.iter_mut().for_each(|matrix| {
        for pos in 0..column {
            if (0..row).any(|cluster| matrix[column * cluster + pos].abs() > 0.001) {
                for cluster in 0..row {
                    matrix[column * cluster + pos] += cv[cluster][pos];
                }
            }
        }
    });
    matrices
}

fn maximize_margin_of(matrices: &[Vec<f64>], row: usize, column: usize) -> Option<(f64, Vec<f64>)> {
    let matrix = matrices
        .iter()
        .map(|matrix| DMatrix::from_row_slice(row, column, &matrix))
        .fold(DMatrix::zeros(column, column), |x, l| {
            let trans = l.transpose();
            let ones = DMatrix::repeat(row, row, 1.);
            let unit = DMatrix::identity(row, row);
            let reg = unit - ones / row as f64;
            x + trans * reg * l
        });
    let eigens = matrix.symmetric_eigen();
    use std::cmp::Ordering::Equal;
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Equal))?;
    if max.1.is_nan() || max.1.is_infinite() {
        None
    } else {
        let argmax = eigens.eigenvectors.column(max.0).iter().copied().collect();
        Some((max.1, argmax))
    }
}

// Return the centrize vector for matrices.
// In other words, it calculates arrays of vector p_1,...,p_{row}, where
// each c_i is column-dimensional vector.
fn centrize_vector_of(matrices: &[Vec<f64>], row: usize, column: usize) -> Vec<Vec<f64>> {
    let mut centrize_vector = vec![vec![0.; column]; row];
    // How many units is considered at position i.
    let mut counts = vec![0; column];
    // Compute the avarage vector at each **column(position)**
    for matrix in matrices.iter() {
        for pos in 0..column {
            if (0..row).any(|r| matrix[column * r + pos].abs() > 0.001) {
                counts[pos] += 1;
                for cluster in 0..row {
                    centrize_vector[cluster][pos] += matrix[column * cluster + pos];
                }
            }
        }
    }
    centrize_vector.iter_mut().for_each(|row_vec| {
        row_vec
            .iter_mut()
            .zip(counts.iter())
            .for_each(|(sum, &count)| *sum /= count as f64)
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

/// Call variant of each pairs of cluster.
/// Note that returned matrix is lower-triangled. In other words,
/// xs[i][j] would be varid only if j < i.
pub fn variant_calling_all_pairs(
    models: &[Vec<POA>],
    data: &[super::Read],
    ws: &[f64],
    c: &super::poa_clustering::ClusteringConfig,
) -> (Vec<Vec<Vec<f64>>>, f64) {
    let (matrices, lk) = calc_matrices_poa(models, data, ws, c);
    let betas: Vec<Vec<Vec<f64>>> = (0..c.cluster_num)
        .map(|i| {
            (0..i)
                .map(|j| call_variants(i, j, &matrices, c.chain_len))
                .collect()
        })
        .collect();
    (betas, lk)
}

pub fn get_lk(
    models: &[Vec<POA>],
    data: &[super::Read],
    c: &super::poa_clustering::ClusteringConfig,
    ws: &[f64],
) -> f64 {
    calc_matrices_poa(models, data, ws, c).1
}

// Call varinants between cluster i and cluster j.
fn call_variants(i: usize, j: usize, matrices: &[Vec<f64>], column: usize) -> Vec<f64> {
    // Extract focal rows for each read.
    let matrices: Vec<Vec<_>> = matrices
        .iter()
        .map(|matrix| {
            let use_pos: Vec<_> = matrix[i * column..(i + 1) * column]
                .iter()
                .zip(matrix[j * column..(j + 1) * column].iter())
                .map(|(&l1, &l2)| l1 > poa_hmm::DEFAULT_LK && l2 > poa_hmm::DEFAULT_LK)
                .collect();
            let class_i = matrix[i * column..(i + 1) * column]
                .iter()
                .zip(use_pos.iter())
                .map(|(&lk, &b)| if b { lk } else { 0. });
            let class_j = matrix[j * column..(j + 1) * column]
                .iter()
                .zip(use_pos.iter())
                .map(|(&lk, &b)| if b { lk } else { 0. });
            class_i.chain(class_j).collect()
        })
        .collect();
    let matrices = centrize(matrices, 2, column);
    match maximize_margin_of(&matrices, 2, column) {
        Some(res) => res.1,
        None => vec![0.; column],
    }
}
