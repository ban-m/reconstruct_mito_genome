extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
#[macro_use]
extern crate log;
extern crate env_logger;
fn main() {
    env_logger::from_env(
        env_logger::Env::default().default_filter_or("variable_chain_length=debug"),
    )
    .init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 3..40;
    let k = 6;
    let len = 150;
    // let coverage = 20;
    let test_num = 100;
    let sample_num: Vec<(usize, u64)> = chain_len
        .flat_map(|e| (0..200).map(|c| (e, c)).collect::<Vec<_>>())
        .collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    println!("HMM\tDist\tCoverage\tLength");
    for coverage in vec![10, 20, 30] {
        let result: Vec<_> = sample_num
            .par_iter()
            .map(|&(chain_len, seed)| {
                debug!("Start{}\t{}", seed, chain_len);
                let (hmm, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
                debug!("Fin{}\t{}", seed, chain_len);
                (hmm, dist, chain_len)
            })
            .collect();
        for (hmm, dist, clen) in result {
            println!("{}\t{}\t{}\t{}", hmm, dist, coverage, len * clen);
        }
    }
}
fn benchmark(
    seed: u64,
    p: &gen_sample::Profile,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    k: usize,
    len: usize,
) -> (f64, u32) {
    let seed = 100342374 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let templates1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates2: Vec<_> = templates1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let prob_1 = 0.5;
    let (_dataset2, model2) = generate_dataset(&templates2, coverage, &mut rng, k);
    let (_dataset1, model1) = generate_dataset(&templates1, coverage, &mut rng, k);
    let dist = templates1
        .iter()
        .zip(templates2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .sum::<u32>();
    let p = &gen_sample::PROFILE;
    let tests: Vec<(_, Vec<_>)> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(prob_1) {
                let d = templates1
                    .iter()
                    .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
                    .collect();
                (1, d)
            } else {
                let d = templates2
                    .iter()
                    .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
                    .collect();
                (2, d)
            }
        })
        .collect();
    debug!("Test...{}\t{}", seed, chain_len);
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = predict(&model1, test);
            let l2 = predict(&model2, test);
            let (l1, l2) = merge_predict(&l1, &l2);
            let p = if l2 < l1 { 1 } else { 2 };
            *ans == p
        })
        .count();
    let hmm = correct as f64 / test_num as f64;
    (hmm, dist)
}

fn generate_dataset<T: Rng>(
    templates: &[Vec<u8>],
    coverage: usize,
    rng: &mut T,
    k: usize,
) -> (Vec<Vec<u8>>, Vec<DBGHMM>) {
    let dataset: Vec<_> = templates
        .iter()
        .map(|e| {
            (0..coverage)
                .map(|_| gen_sample::introduce_randomness(e, rng, &gen_sample::PROFILE))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut f = Factory::new();
    let models: Vec<_> = dataset
        .iter()
        .map(|e| {
            let ds:Vec<_> = e.iter().map(|e|e.as_slice()).collect();
            let ws = vec![1.; e.len()];
            f.generate_with_weight(&ds, &ws, k,&mut vec![])
        })
        .collect();
    let dataset: Vec<_> = (0..coverage)
        .map(|i| {
            dataset
                .iter()
                .flat_map(|e| &e[i])
                .copied()
                .collect::<Vec<_>>()
        })
        .collect();
    (dataset, models)
}

fn predict(models: &[DBGHMM], test: &[Vec<u8>]) -> Vec<f64> {
    models
        .iter()
        .zip(test.iter())
        // .map(|(e, f)| e.forward(f, &PACBIO_CONFIG))
        .map(|(e, f)| e.forward(f, &DEFAULT_CONFIG))
        .collect()
}

fn merge_predict(l1: &[f64], l2: &[f64]) -> (f64, f64) {
    // l1 and l2 are the log prob.
    let sum1 = l1.iter().sum();
    let sum2 = l2.iter().sum();
    (sum1, sum2)
}
