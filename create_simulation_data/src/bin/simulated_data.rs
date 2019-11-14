extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
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
            let (hmm, aln, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
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
    let (dataset2, model2) = generate_dataset(&templates2, coverage, &mut rng, k);
    let (dataset1, model1) = generate_dataset(&templates1, coverage, &mut rng, k);
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
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let test = test.concat();
            let l1: u32 = dataset1
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, &test))
                .min()
                .unwrap();
            let l2: u32 = dataset2
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, &test))
                .min()
                .unwrap();
            if l2 < l1 {
                *ans == 2
            } else if l1 < l2 {
                *ans == 1
            } else {
                ans == [1, 2].choose(&mut rng).unwrap()
            }
        })
        .count();
    let aln = correct as f64 / test_num as f64;
    (hmm, aln, dist)
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
        .collect()
}

#[allow(dead_code)]
fn merge_predict_naive(l1: &[f64], l2: &[f64]) -> (f64, f64) {
    // l1 and l2 are the log prob.
    let sum1 = l1.iter().sum();
    let sum2 = l2.iter().sum();
    as_weight(sum1, sum2)
}

#[allow(dead_code)]
fn merge_predict(l1: &[f64], l2: &[f64]) -> (f64, f64) {
    //println!("merging prediction below:");
    let ratio: Vec<_> = l1
        .iter()
        .zip(l2.iter())
        .map(|(&x1, &x2)| as_weight(x1, x2))
        .collect();
    let weights: Vec<_> = ratio
        .iter()
        .map(|(x1, x2)| 2f64.ln() + x1 * x1.ln() + x2 * x2.ln())
        .collect();
    let tot = weights.iter().fold(0., |x, y| x + y);
    // eprintln!("Dump weights");
    // for (idx, ((&l1, &l2), w)) in l1.iter().zip(l2.iter()).zip(weights.iter()).enumerate() {
    //     let (w1, w2) = as_weight(l1, l2);
    //     eprintln!("{}\t{:.4}\t{:.4}\t{:.4}", idx, w1, w2, w / tot);
    // }
    let (p1, p2) = ratio
        .into_iter()
        .zip(weights.into_iter())
        .map(|((f1, f2), w)| (f1 * w, f2 * w))
        .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
    assert!((p1 / tot + p2 / tot - 1.).abs() < 0.0001);
    // eprintln!("Summping to:{}\t{}",p1 / tot, p2 / tot);
    // eprintln!();
    (p1 / tot, p2 / tot)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}
