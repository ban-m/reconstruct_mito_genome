extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
fn main() {
    use env_logger::Env;
    env_logger::Builder::from_env(Env::default().default_filter_or("em_algorithm_check=debug"))
        .init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let coverage = 6;
    let test_num = 150;
    let sample_num: Vec<_> = (0..1).collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    debug!("Start");
    let result: Vec<_> = sample_num
        .into_iter()
        .map(|seed| {
            let (em, naive, testnum) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
            (em, naive, testnum, coverage)
        })
        .collect();
    println!("EM\tNaive\tCoverage\tTestNum\tLength");
    for (em, naive, tn, coverage) in result {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            em,
            naive,
            coverage,
            tn,
            len * chain_len
        );
        debug!(
            "{:.3}\t{:.3}",
            em as f64 / tn as f64,
            naive as f64 / tn as f64
        )
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
) -> (usize, usize, usize) {
    let seed = 100342374 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let templates1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates2: Vec<_> = templates1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    debug!("Generateing dataset...");
    let (dataset, assignment, answer, border) =
        generate_dataset(&templates1, &templates2, coverage, test_num, &mut rng);
    debug!("Gen Dataset");
    let naive_pred = naive_solve(&dataset, &assignment, border, k);
    let em_pred = em_solve(&dataset, &assignment, border, k, seed);
    assert_eq!(em_pred.len(), answer.len());
    let em = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|&(&p, &a)| p == a)
        .count();
    let naive = naive_pred
        .iter()
        .zip(answer.iter())
        .filter(|&(&p, &a)| p == a)
        .count();
    (em, naive, answer.len())
}

fn generate_dataset<T: Rng>(
    template1: &[Vec<u8>],
    template2: &[Vec<u8>],
    coverage: usize,
    test_num: usize,
    rng: &mut T,
) -> (Vec<Vec<Vec<u8>>>, Vec<u8>, Vec<u8>, usize) {
    let answer: Vec<_> = (0..test_num + 2 * coverage)
        .map(|i| if i < 4 { i % 2 == 0 } else { rng.gen_bool(0.2) })
        .map(|b| if b { 0 } else { 1 })
        .collect();
    let mut gen = |t: &[Vec<u8>]| {
        t.iter()
            .map(|e| gen_sample::introduce_randomness(e, rng, &gen_sample::PROFILE))
            .collect::<Vec<_>>()
    };
    debug!("Coverage:{}\ttest_num:{}", coverage, test_num);
    let dataset: Vec<_> = answer
        .iter()
        .map(|e| {
            if e % 2 == 0 {
                gen(template1)
            } else {
                gen(template2)
            }
        })
        .collect();
    let assignment: Vec<_> = answer.iter().take(2 * coverage).copied().collect();
    let answer: Vec<_> = answer.into_iter().skip(2 * coverage).collect();
    assert_eq!(dataset.len(), assignment.len() + answer.len());
    let l = assignment.len();
    (dataset, assignment, answer, l)
}

fn naive_solve(dataset: &[Vec<Vec<u8>>], label: &[u8], border: usize, k: usize) -> Vec<u8> {
    let (d1, d2): (Vec<_>, Vec<_>) = dataset
        .iter()
        .take(border)
        .zip(label)
        .partition(|&(_, &ans)| ans == 0);
    let d1: Vec<_> = d1.into_iter().map(|(x, _)| x).collect();
    let d2: Vec<_> = d2.into_iter().map(|(x, _)| x).collect();
    let model1 = construct_from_reads(d1, k);
    let model2 = construct_from_reads(d2, k);
    dataset
        .par_iter()
        .skip(border)
        .map(|read| {
            let lkdiff = read
                .iter()
                .enumerate()
                .map(|(idx, chunk)| {
                    model1[idx].forward(chunk, &DEFAULT_CONFIG)
                        - model2[idx].forward(chunk, &DEFAULT_CONFIG)
                })
                .sum::<f64>();
            // debug!("PRED:{}", lkdiff);
            if lkdiff.is_sign_positive() {
                0
            } else {
                1
            }
        })
        .collect()
}

fn em_solve(dataset: &[Vec<Vec<u8>>], label: &[u8], border: usize, k: usize, s: u64) -> Vec<u8> {
    let pred = naive_solve(dataset, label, border, k);
    debug!("Bootstrapped by Niave approach");
    let (mut w0, mut w1) = {
        let (w0, w1) = label.iter().chain(pred.iter()).fold((0, 0), |(x, y), &p| {
            if p == 0 {
                (x + 1, y)
            } else {
                (x, y + 1)
            }
        });
        let len = dataset.len() as f64;
        (w0 as f64 / len, w1 as f64 / len)
    };
    debug!("(w0,w1)=({},{})", w0, w1);
    let mut gamma0: Vec<_> = label
        .iter()
        .chain(pred.iter())
        .map(|&p| if p == 0 { 1. } else { 0. })
        .collect();
    let mut gamma1: Vec<_> = label
        .iter()
        .chain(pred.iter())
        .map(|&p| if p == 0 { 0. } else { 1. })
        .collect();
    debug!("Calculated gammas.");
    let mut model0 = construct_with_weights(dataset, &gamma0, k);
    let mut model1 = construct_with_weights(dataset, &gamma1, k);
    let mut lk = calc_lk(&dataset[border..], &model0, &model1, w0, w1);
    debug!("Current log likelihood:{:.3}", lk);
    //    use std::io::{BufWriter, Write};
    let ep = 0.0001;
    for _ in 0..10 {
        for (idx, read) in dataset.iter().enumerate().skip(border) {
            let log_m0 = read
                .iter()
                .enumerate()
                .map(|(idx, chunk)| model0[idx].forward(chunk, &DEFAULT_CONFIG))
                .sum::<f64>();
            let log_m1 = read
                .iter()
                .enumerate()
                .map(|(idx, chunk)| model1[idx].forward(chunk, &DEFAULT_CONFIG))
                .sum::<f64>();
            let w = calc_logsum(log_m0, log_m1, w0, w1);
            gamma0[idx] = (w0.ln() + log_m0 - w).exp();
            gamma1[idx] = (w1.ln() + log_m1 - w).exp();
        }
        let tot = gamma0.iter().sum::<f64>() + gamma1.iter().sum::<f64>();
        w0 = gamma0.iter().sum::<f64>() / tot;
        w1 = gamma1.iter().sum::<f64>() / tot;
        model0 = construct_with_weights(&dataset, &gamma0, k);
        model1 = construct_with_weights(&dataset, &gamma1, k);
        let next_lk = calc_lk(&dataset[border..], &model0, &model1, w0, w1);
        debug!("Updated LK:{:.3}", next_lk);
        debug!("(w0,w1)=({},{})", w0, w1);
        if (next_lk - lk).abs() < ep {
            break;
        } else {
            lk = next_lk;
        }
    }
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(s + 9893242349);
    gamma0
        .iter()
        .skip(border)
        .map(|&g0| if rng.gen_range(0., 1.) < g0 { 0 } else { 1 })
        .collect()
}

fn calc_lk(ds: &[Vec<Vec<u8>>], m0: &[DBGHMM], m1: &[DBGHMM], w0: f64, w1: f64) -> f64 {
    ds.iter()
        .map(|read| {
            read.par_iter()
                .enumerate()
                .map(|(idx, chunk)| {
                    let m0 = m0[idx].forward(chunk, &DEFAULT_CONFIG);
                    let m1 = m1[idx].forward(chunk, &DEFAULT_CONFIG);
                    calc_logsum(m0, m1, w0, w1)
                })
                .sum::<f64>()
        })
        .sum::<f64>()
}

fn calc_logsum(l0: f64, l1: f64, w0: f64, w1: f64) -> f64 {
    // Log (w0 * exp(l0) + w1 * exp(l1))
    let f1 = w0.ln() + l0;
    let f2 = w1.ln() + l1;
    let m = f1.max(f2);
    m + ((f1 - m).exp() + (f2 - m).exp()).ln()
}

fn construct_from_reads(ds: Vec<&Vec<Vec<u8>>>, k: usize) -> Vec<DBGHMM> {
    let len = ds[0].len();
    let mut res: Vec<Vec<&[u8]>> = vec![vec![]; len];
    assert!(ds.iter().all(|e| e.len() == len));
    for read in ds {
        for (idx, chunk) in read.iter().enumerate() {
            res[idx].push(chunk);
        }
    }
    let mut f = Factory::new();
    res.into_iter()
        .map(|e| f.generate_from_ref(&e, k))
        .collect()
}

fn construct_with_weights(ds: &[Vec<Vec<u8>>], ws: &[f64], k: usize) -> Vec<DBGHMM> {
    let len = ds[0].len();

    let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; len];
    for read in ds.into_iter() {
        for (idx, chunk) in read.iter().enumerate() {
            chunks[idx].push(chunk);
        }
    }
    let mut f = Factory::new();
    chunks
        .into_iter()
        .map(|cs| f.generate_with_weight(&cs, &ws, k))
        .collect()
}
