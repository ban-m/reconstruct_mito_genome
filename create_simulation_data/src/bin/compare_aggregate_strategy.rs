extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const BADREAD_CONFIG: dbg_hmm::Config = dbg_hmm::Config {
    mismatch: 0.0344,
    p_match: 0.88,
    p_ins: 0.0549,
    p_del: 0.0651,
    p_extend_ins: 0.0337,
    p_extend_del: 0.1787,
    p_del_to_ins: 0.0,
    match_score: 1,
    mism_score: -1,
    del_score: -1,
    ins_score: -1,
    base_freq: [0.2543, 0.2290, 0.2558, 0.2609],
};
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let k = 6;
    let len = 100;
    let min_coverage = 10;
    let max_coverage = 31;
    let test_num = 100;
    let sample_num: Vec<(u64, usize)> = (0..500)
        .flat_map(|e| {
            (min_coverage..max_coverage)
                .map(|c| (e, c))
                .collect::<Vec<_>>()
        })
        .collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let result: Vec<_> = sample_num
        .into_par_iter()
        .map(|(seed, coverage)| {
            let (na, w, wl, ind, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
            (na, w, wl, ind, dist, coverage)
        })
        .collect();
    println!("Niave\tWeightedSum\tWeightedLikelihood\tIndependent\tDist\tCoverage\tLength");
    for (w, x, y, z, dist, coverage) in result {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            w,
            x,
            y,
            z,
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
) -> (f64, f64, f64, f64, u32) {
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
    let dist = templates1
        .iter()
        .zip(templates2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .sum::<u32>();
    let (_dataset2, model2) = generate_dataset(&templates2, coverage, &mut rng, k);
    let (_dataset1, model1) = generate_dataset(&templates1, coverage, &mut rng, k);
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
    let add = |x: f64, a: bool| if a { x + 1. } else { x };
    let (w, x, y, z) = tests
        .iter()
        .map(|&(ans, ref test)| {
            let l1 = predict(&model1, test);
            let l2 = predict(&model2, test);
            // x,y, and z are 1 or 2.
            let (w, x, y, z) = merge_predict(&l1, &l2);
            (ans == w, ans == x, ans == y, ans == z)
        })
        .fold((0., 0., 0., 0.), |(w, x, y, z), (a, b, c, d)| {
            (add(w, a), add(x, b), add(y, c), add(z, d))
        });
    let d = test_num as f64;
    let (w, x, y, z) = (w / d, x / d, y / d, z / d);
    (w, x, y, z, dist)
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
    let models: Vec<_> = dataset.iter().map(|e| f.generate(e, k)).collect();
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
        // .map(|(e, f)| e.forward(f, &BADREAD_CONFIG))
        .collect()
}

fn merge_predict(l1: &[f64], l2: &[f64]) -> (u8, u8, u8, u8) {
    // naive calc. In other words, we just multiple all the likelihood.
    // Note that l1 and l2 are log-likelihood thus we just add all the values.
    let naive = l1
        .iter()
        .zip(l2.iter())
        .map(|(l1, l2)| l1 - l2)
        .sum::<f64>();
    let naive = if naive.is_sign_positive() { 1 } else { 2 };
    // Weighted sum. First convert l1,l2 to the ratios.
    let ratio: Vec<_> = l1
        .iter()
        .zip(l2.iter())
        .map(|(&x1, &x2)| as_weight(x1, x2))
        .collect();
    let weights: Vec<_> = ratio
        .iter()
        .map(|(x1, x2)| 2f64.ln() + x1 * x1.ln() + x2 * x2.ln())
        .collect();
    let s_o_w = ratio
        .iter()
        .map(|(w1, w2)| w1 - w2)
        .zip(weights.iter())
        .map(|(x, y)| x * y)
        .sum::<f64>();
    let s_o_w = if s_o_w.is_sign_positive() { 1 } else { 2 };
    // Weighted log likelihood.
    let s_o_l = l1
        .iter()
        .zip(l2.iter())
        .map(|(l1, l2)| l1 - l2)
        .zip(weights.iter())
        .map(|(x, y)| x * y)
        .sum::<f64>();
    let s_o_l = if s_o_l.is_sign_positive() { 1 } else { 2 };
    // Independent prediction.
    let ind = l1
        .iter()
        .zip(l2.iter())
        .map(|(l1, l2)| if l1 > l2 { 1 } else { -1 })
        .sum::<i32>();
    let ind = if ind.is_positive() { 1 } else { 2 };
    (naive, s_o_w, s_o_l, ind)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}
