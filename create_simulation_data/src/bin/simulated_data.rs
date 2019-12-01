extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use create_simulation_data::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let min_coverage = 6;
    let max_coverage = 40;
    let test_num = 100;
    let sample_num: Vec<(u64, usize)> = (0..100)
        .flat_map(|e| {
            (min_coverage..max_coverage)
                .step_by(4)
                .map(|c| (e, c))
                .collect::<Vec<_>>()
        })
        .collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    use std::time::Instant;
    let result: Vec<_> = sample_num
        .into_par_iter()
        .map(|(seed, coverage)| {
            let s = Instant::now();
            let (hmm, aln, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
            eprintln!("{:?}", Instant::now() - s);
            (hmm, aln, dist, coverage)
        })
        .collect();
    println!("HMM\tAln\tDist\tCoverage\tLength");
    for (hmm, aln, dist, coverage) in result {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            hmm,
            aln,
            dist,
            coverage,
            len * chain_len
        );
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
) -> (f64, f64, u32) {
    let seed = 100342374 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let template1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let template2: Vec<_> = template1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let dist = template1
        .iter()
        .zip(template2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .sum::<u32>();
    let prob_0 = 0.5;
    let (dataset, label, answer, border) =
        generate_dataset(&template1, &template2, coverage, test_num, &mut rng, prob_0);
    let em_pred = em_solve(&dataset, &label, border, k, &answer);
    let correct = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|(p, a)| a == p)
        .count();
    let hmm = correct as f64 / test_num as f64;
    let aln = align_solve(&dataset, &label, border);
    let correct = aln
        .iter()
        .zip(answer.iter())
        .filter(|(p, a)| p == a)
        .count();
    let aln = correct as f64 / test_num as f64;
    (hmm, aln, dist)
}
