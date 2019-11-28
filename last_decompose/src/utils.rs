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
        .fold(f64::MIN, |x, y| x.min(y));
    m + ls
        .iter()
        .zip(ws.iter())
        .map(|(&l, &w)| l + w.ln() - m)
        .map(|e| e.exp())
        .sum::<f64>()
        .ln()
}

