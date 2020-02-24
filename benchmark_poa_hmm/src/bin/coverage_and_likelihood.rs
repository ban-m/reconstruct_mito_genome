extern crate edlib_sys;
extern crate env_logger;
extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::DEFAULT_CONFIG;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 150;
    println!("Seed\tCoverage\tLikelihoodRatio\tOriginalLK\tNumEdges");
    let rep = 20;
    let covs: Vec<_> = (15..num_seq).collect();
    let reps: Vec<_> = (0..rep).collect();
    use poa_hmm::gen_sample::*;
    let result: Vec<_> = reps
        .into_par_iter()
        .flat_map(|seed| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(121_892 + seed);
            let template: Vec<_> = poa_hmm::gen_sample::generate_seq(&mut rng, len);
            let data: Vec<Vec<_>> = (0..num_seq)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let max = {
                let cov = 500;
                let data: Vec<Vec<_>> = (0..cov)
                    .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                    .collect();
                let m = poa_hmm::POA::generate_vec(&data);
                data.iter()
                    .map(|q| m.forward(&q, &DEFAULT_CONFIG))
                    .sum::<f64>()
                    / cov as f64
            };
            let mut result = vec![];
            for &cov in &covs {
                for chunk in data.chunks_exact(cov) {
                    let m = poa_hmm::POA::generate_vec(chunk);
                    for q in chunk {
                        let lk = m.forward(&q, &DEFAULT_CONFIG);
                        result.push((seed, cov, max - lk, lk))
                    }
                }
            }
            result
        })
        .collect();
    for (seed, cov, ratio, orig) in result {
        println!("{}\t{}\t{}\t{}", seed, cov, ratio, orig);
    }
}
