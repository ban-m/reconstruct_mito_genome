extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 150;
    let k = 6;
    println!("Seed\tCoverage\tLikelihoodRatio\tOriginalLK\tNumEdges");
    let rep = 20;
    let covs: Vec<_> = (15..num_seq).collect();
    let reps: Vec<_> = (0..rep).collect();
    let result: Vec<_> = reps
        .into_par_iter()
        .flat_map(|seed| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(121_892 + seed);
            let template: Vec<_> = generate_seq(&mut rng, len);
            let data: Vec<Vec<_>> = (0..num_seq)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let mut f = Factory::new();
            let max = {
                let cov = 500;
                let data: Vec<Vec<_>> = (0..cov)
                    .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                    .collect();
                let data: Vec<_> = data.iter().map(|e| e.as_slice()).collect();
                let m = f.generate_with_weight(&data, &vec![1.; cov], k, &mut vec![]);
                data.iter()
                    .map(|q| m.forward(&q, &DEFAULT_CONFIG))
                    .sum::<f64>()
                    / cov as f64
            };
            let mut result = vec![];
            for &cov in &covs {
                let w = vec![1.; cov];
                for chunk in data.chunks_exact(cov) {
                    let m: Vec<_> = chunk.iter().map(|e| e.as_slice()).collect();
                    let m = f.generate_with_weight(&m, &w, k, &mut vec![]);
                    let edges = m.edge_num();
                    for q in chunk {
                        let lk = m.forward(&q, &DEFAULT_CONFIG);
                        result.push((seed, cov, max - lk, lk, edges))
                    }
                }
            }
            result
        })
        .collect();
    for (seed, cov, ratio, orig, edge) in result {
        println!("{}\t{}\t{}\t{}\t{}", seed, cov, ratio, orig, edge);
    }
}
