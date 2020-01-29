use dbg_hmm::Config;
use dbg_hmm::DBGHMM;
use na::{DMatrix, DVector};
// Retun a vector with the length of 1 s.t. the margin would be maximized.
fn variant_call(models: &[Vec<Vec<DBGHMM>>], data: &[ERead], ws: &[f64], c: &Config) -> Vec<f64> {
    let num_cluster = models.len();
    let num_chain = models[0].iter().map(|e|e.len()).sum:<usize>();
    let reads: Vec<_> = data.iter()
        .flat_map(|i| models.iter().flat_map(|ms|
        ))
        .map(|data| DMatrix::from_column_slice(num_cluster, num_chain, &data))
        .collect();
    let data = reads
        .into_iter()
        .fold(DMatrix::zeros(num_chain, num_chain), |x, l| {
            let trans = l.clone().transpose();
            let ltl = trans.clone() * l.clone();
            let ones = DMatrix::repeat(num_cluster, num_cluster, 1.);
            let lt_one_l = trans * ones * l;
            x + ltl - lt_one_l / num_cluster as f64
        });
    let eigens = data.clone().symmetric_eigen();
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    eprintln!("{:?}", max);
    let eigenvec = eigens.eigenvectors.column(max.0);
    let answer = DVector::from_column_slice(&pos);
    eprintln!("{}\n{}", eigenvec, answer);
    // eprintln!("{}", data * eigenvec);
    // eprintln!("{}", max.1 * eigenvec);
}
