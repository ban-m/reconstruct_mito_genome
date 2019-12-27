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
    println!("Seed\tCoverage\tContainedSelf\tNotContained");
    let rep = 15;
    let covs: Vec<_> = (1..num_seq).collect();
    let result = (0..rep).collect::<Vec<_>>();
    let result: Vec<_> = result
        .iter()
        .flat_map(|&seed| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(seed as u64);
            let template: Vec<_> = generate_seq(&mut rng, len);
            let data: Vec<Vec<_>> = (0..num_seq)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let tests: Vec<_> = (0..num_seq)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let mut f = Factory::new();
            let weight = vec![1.; num_seq];
            covs.iter()
                .flat_map(|&coverage| {
                    data.chunks_exact(coverage)
                        .zip(tests.chunks_exact(coverage))
                        .flat_map(|(data, test)| {
                            let m1: Vec<_> = data.iter().map(|e| e.as_slice()).collect();
                            let m1 = f.generate_with_weight_prior(
                                &m1,
                                &weight[..coverage],
                                k,
                                &mut vec![],
                            );
                            let m2: Vec<_> = test.iter().map(|e| e.as_slice()).collect();
                            let m2 = f.generate_with_weight_prior(
                                &m2,
                                &weight[..coverage],
                                k,
                                &mut vec![],
                            );
                            test.iter()
                                .map(|t| {
                                    let m1 = m1.forward(t, &DEFAULT_CONFIG);
                                    let m2 = m2.forward(t, &DEFAULT_CONFIG);
                                    (seed, coverage, m2, m1)
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect();
    for (seed, cov, contained, not_contained) in result {
        println!("{}\t{}\t{}\t{}", seed, cov, contained, not_contained);
    }
}
