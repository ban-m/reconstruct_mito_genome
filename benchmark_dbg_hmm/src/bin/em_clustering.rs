extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let len = 150;
    let k = 6;
    let num_seq = (0..10).collect::<Vec<usize>>();
    let test_num = (10..400).step_by(15).collect::<Vec<usize>>();
    let sample_num: Vec<u64> = (0..50).collect();
    let p = &gen_sample::Profile {
        sub: 0.003,
        ins: 0.004,
        del: 0.003,
    };
    let c = if args[1] == "default" {
        DEFAULT_CONFIG
    } else {
        PACBIO_CONFIG
    };
    let s = &gen_sample::PROFILE;
    let result: Vec<_> = sample_num
        .par_iter()
        .flat_map(|&seed| {
            let mut res = vec![];
            for training in &num_seq {
                for test in &test_num {
                    let (n, e) = benchmark(p, s, seed, k, len, *training, *test, &c);
                    res.push((training, test, n, e));
                }
            }
            eprintln!("Fin{}", seed);
            res
        })
        .collect();
    println!("Training\tTes\tNaive\tEM");
    for (training, test, n, e) in result {
        println!("{}\t{}\t{}\t{}", training, test, n, e);
    }
}
fn benchmark(
    p: &gen_sample::Profile,
    s: &gen_sample::Profile,
    seed: u64,
    k: usize,
    len: usize,
    training: usize,
    test_num: usize,
    config: &Config,
) -> (f64, f64) {
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(122492 + seed);
    let template1 = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let template2 = dbg_hmm::gen_sample::introduce_randomness(&template1, &mut rng, &p);
    let answers: Vec<_> = (0..(training + test_num))
        .map(|_| rng.gen_bool(0.5))
        .collect();
    let dataset: Vec<Vec<_>> = answers
        .iter()
        .map(|&e| {
            if e {
                gen_sample::introduce_randomness(&template1, &mut rng, s)
            } else {
                gen_sample::introduce_randomness(&template2, &mut rng, s)
            }
        })
        .collect();
    let pred = naive_pred(&dataset, &answers[..training], k, config);
    let (mut w0, mut w1) = {
        let w0 = answers.iter().filter(|&&e| e).count();
        let w1 = answers.iter().filter(|&&e| !e).count();
        let len = answers.len() as f64;
        (w0 as f64 / len, w1 as f64 / len)
    };
    let mut gamma0: Vec<_> = pred.iter().map(|&e| if e { 1. } else { 0. }).collect();
    let mut gamma1: Vec<_> = pred.iter().map(|&e| if !e { 1. } else { 0. }).collect();
    let mut model0 = construct_model(&dataset, &gamma0, k);
    let mut model1 = construct_model(&dataset, &gamma1, k);
    let mut lks = vec![];
    let mut lk = calc_lk(&dataset[training..], &model0, &model1, w0, w1, config);
    let ep = 0.001;
    for _ in 0..10 {
        lks.push(lk);
        for (idx, chunk) in dataset.iter().enumerate().skip(training) {
            let log_m0 = model0.forward(chunk, config);
            let log_m1 = model1.forward(chunk, config);
            let w = logsum(log_m0, log_m1, w0, w1);
            gamma0[idx] = (w0.ln() + log_m0 - w).exp();
            gamma1[idx] = (w1.ln() + log_m1 - w).exp();
        }
        let tot = gamma0.iter().chain(gamma1.iter()).sum::<f64>();
        w0 = gamma0.iter().sum::<f64>() / tot;
        w1 = gamma1.iter().sum::<f64>() / tot;
        model0 = construct_model(&dataset, &gamma0, k);
        model1 = construct_model(&dataset, &gamma0, k);
        let next_lk = calc_lk(&dataset[training..], &model0, &model1, w0, w1, config);
        if (lk - next_lk).abs() < ep {
            break;
        } else {
            lk = next_lk;
        }
    }
    let em_pred: Vec<_> = gamma0
        .iter()
        .zip(gamma1.iter())
        .map(|(f1, f2)| f1 > f2)
        .collect();
    let len = answers.len() as f64;
    let em = em_pred
        .iter()
        .zip(answers.iter())
        .filter(|(e, f)| e == f)
        .count();
    let naive = pred
        .iter()
        .zip(answers.iter())
        .filter(|(e, f)| e == f)
        .count();
    (em as f64 / len, naive as f64 / len)
}

fn naive_pred(dataset: &[Vec<u8>], label: &[bool], k: usize, c: &Config) -> Vec<bool> {
    if label.is_empty() {
        dataset
            .iter()
            .enumerate()
            .map(|(idx, _)| idx % 2 == 0)
            .collect()
    } else {
        let training = label.len();
        let m0: Vec<_> = dataset
            .iter()
            .zip(label.iter())
            .filter(|&(_, &b)| b)
            .map(|(e, _)| e.as_slice())
            .collect();
        let m1: Vec<_> = dataset
            .iter()
            .zip(label.iter())
            .filter(|&(_, &b)| !b)
            .map(|(e, _)| e.as_slice())
            .collect();
        let m0 = DBGHMM::new_from_ref(&m0, k);
        let m1 = DBGHMM::new_from_ref(&m1, k);
        let pred = dataset
            .iter()
            .skip(training)
            .map(|chunk| m0.forward(chunk, c) > m1.forward(chunk, c));
        label.iter().copied().chain(pred).collect()
    }
}

fn construct_model(dataset: &[Vec<u8>], ws: &[f64], k: usize) -> DBGHMM {
    let mut f = Factory::new();
    let d: Vec<_> = dataset.iter().map(|e| e.as_slice()).collect();
    f.generate_with_weight(&d, ws, k)
}

fn calc_lk(dataset: &[Vec<u8>], m0: &DBGHMM, m1: &DBGHMM, w0: f64, w1: f64, c: &Config) -> f64 {
    dataset
        .iter()
        .map(|chunk| {
            let m0 = m0.forward(chunk, c);
            let m1 = m1.forward(chunk, c);
            logsum(m0, m1, w0, w1)
        })
        .sum::<f64>()
}

fn logsum(l0: f64, l1: f64, w0: f64, w1: f64) -> f64 {
    // Log (w0 * exp(l0) + w1 * exp(l1))
    let f1 = w0.ln() + l0;
    let f2 = w1.ln() + l1;
    let m = f1.max(f2);
    m + ((f1 - m).exp() + (f2 - m).exp()).ln()
}
