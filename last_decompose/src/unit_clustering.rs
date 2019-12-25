//! Unit level clustering(experimental)
use super::*;

fn initial_weights(datanum: usize, cluster_num: usize) -> Vec<Vec<f64>> {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(datanum as u64 * 21);
    let num_of_ball = cluster_num * NUM_OF_BALL;
    let denom = (num_of_ball as f64).recip();
    let gen_dist = |_| {
        let mut bucket = vec![0; cluster_num];
        (0..num_of_ball).for_each(|_| bucket[rng.gen_range(0, cluster_num)] += 1);
        bucket.iter().map(|&e| e as f64 * denom).collect::<Vec<_>>()
    };
    let weights: Vec<Vec<_>> = (0..datanum).map(gen_dist).collect();
    weights
}

pub fn unit_clustering(data: &[&[u8]], k: usize, cluster_num: usize, c: &Config) -> Vec<u8> {
    let (mut f, mut buf) = (Factory::new(), vec![]);
    let mut weights = initial_weights(data.len(), cluster_num);
    let mut rhos = weights.clone();
    let mut models: Vec<_> = (0..cluster_num)
        .map(|cl| {
            let ws: Vec<_> = weights.iter().map(|e| e[cl]).collect();
            f.generate_with_weight_prior(data, &ws, k, &mut buf)
        })
        .collect();
    let mut alphas: Vec<_> = (0..cluster_num)
        .map(|cl| weights.iter().map(|ws| ws[cl]).sum::<f64>())
        .collect();
    // 350 loops
    let betas = (0..)
        .map(|i| 0.01 * 1.02_f64.powi(i as i32))
        .take_while(|&e| e <= 1.)
        .chain(vec![1.]);
    let mut cluster_weight = vec![vec![0.; data.len()]; cluster_num];
    for beta in betas {
        let mut soe = 100.;
        let mut soe_diff = 100.;
        for _ in 0..100 {
            if soe < 0.011 || soe_diff < 0.11 {
                break;
            }
            let tot = alphas.iter().sum::<f64>();
            data.par_iter()
                .zip(weights.par_iter_mut())
                .zip(rhos.par_iter_mut())
                .for_each(|((data, weights), log_rhos)| {
                    models
                        .iter()
                        .zip(alphas.iter())
                        .zip(log_rhos.iter_mut())
                        .for_each(|((model, a), r)| {
                            let weight = model.weight();
                            let coverage_offset = offset(weight, A_PRIOR, B_PRIOR);
                            let lk = model.forward(data, c);
                            let lk = lk + coverage_offset;
                            let lk = lk.max(c.null_model(data));
                            *r = beta * (lk + digamma(*a) - digamma(tot))
                        });
                    let log_sum_rho = utils::logsumexp(&log_rhos);
                    weights
                        .iter_mut()
                        .zip(log_rhos)
                        .for_each(|(w, r)| *w = (*r - log_sum_rho).exp());
                });
            models
                .iter_mut()
                .enumerate()
                .zip(cluster_weight.iter_mut())
                .for_each(|((cl, m), ws)| {
                    ws.iter_mut()
                        .zip(weights.iter())
                        .for_each(|(w, r)| *w = r[cl]);
                    *m = f.generate_with_weight_prior(data, &ws, k, &mut buf)
                });
            alphas = (0..cluster_num)
                .map(|cl| weights.iter().map(|ws| ws[cl]).sum::<f64>() + ALPHA)
                .map(|alpha| (alpha - 1.) * beta + 1.)
                .collect();
            let next_soe = weights.iter().map(|e| entropy(e)).sum::<f64>();
            soe_diff = soe - next_soe;
            soe = next_soe;
        }
        report(soe, &cluster_weight, beta);
    }
    weights
        .iter()
        .map(|ws| {
            ws.iter()
                .enumerate()
                .max_by(|e, f| (e.1).partial_cmp(f.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx as u8)
                .unwrap_or(0)
        })
        .collect()
}
fn report(soe: f64, cw: &[Vec<f64>], beta: f64) {
    let pi: Vec<_> = cw
        .iter()
        .map(|e| e.iter().sum::<f64>())
        .map(|e| format!("{:.2}", e))
        .collect();
    let pi = pi.join("\t");
    info!("Summary\t{:.3}\t{:.3}\t{}", beta, soe, pi,);
}
