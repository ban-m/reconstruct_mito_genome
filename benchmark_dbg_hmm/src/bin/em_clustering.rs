extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
use std::time::Instant;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let len = 150;
    let k = 6;
    // (num_seq, test_num, sample_num);
    let params: Vec<_> = (10..11)
        .flat_map(|x| {
            (20..21)
                .step_by(15)
                .flat_map(|y| (0..20).map(|z| (x, y, z)).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        })
        .collect();
    let p = &gen_sample::Profile {
        sub: 0.006,
        ins: 0.006,
        del: 0.006,
    };
    let c = if args[1] == "default" {
        DEFAULT_CONFIG
    } else {
        PACBIO_CONFIG
    };
    let s = &gen_sample::PROFILE;
    println!("Training\tTest\tNaive\tEM\tDist");
    let res: Vec<_> = params
        .into_par_iter()
        .map(|(training, test, seed)| {
            let start = Instant::now();
            info!("Start:{}\t{}\t{}", training, test, seed);
            let (e, n, d) = benchmark(p, s, seed, k, len, training, test, &c);
            let ela = Instant::now() - start;
            let millis = ela.as_millis() as f64 / (training + test) as f64;
            info!(
                "End {}\t{}\t{}, {:?}({:.4}/sample)",
                training, test, seed, ela, millis
            );
            (training, test, n, e, d)
        })
        .collect();
    for (training, test, n, e, d) in res {
        println!("{}\t{}\t{:.4}\t{:.4}\t{}", training, test, n, e, d);
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
) -> (f64, f64, u32) {
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(138222324 + seed);
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
    let d = edlib_sys::global_dist(&template1, &template2);
    let pred = naive_pred(&dataset, &answers, training, k, config);
    // assert!(dataset.len() == pred.len());
    let em_pred = em_pred(&dataset, &answers[..training], k, config);
    let naive = pred
        .iter()
        .zip(answers.iter())
        .skip(training)
        .filter(|(e, f)| e == f)
        .count() as f64;
    let em = em_pred
        .iter()
        .zip(answers.iter())
        .skip(training)
        .filter(|(e, f)| e == f)
        .count() as f64;
    (em / test_num as f64, naive / test_num as f64, d)
}

fn em_pred(dataset: &[Vec<u8>], label: &[bool], k: usize, c: &Config) -> Vec<bool> {
    let label: Vec<_> = label.iter().map(|&e| if e { 0 } else { 1 }).collect();
    let data:Vec<_> = dataset
        .iter()
        .map(|chunk| {
            let read = vec![chunk.to_vec()];
            last_decompose::ERead::new_with_lowseq(read, "test")
        })
        .collect();
    let forbidden = vec![vec![]; data.len()];
    let answer = vec![0; data.len()];
    last_decompose::clustering(&data, &label, &forbidden, k, 2, &[1], &answer, c)
        .into_iter()
        .map(|e| e == 0)
        .collect()
}

fn naive_pred(
    dataset: &[Vec<u8>],
    answer: &[bool],
    training: usize,
    _k: usize,
    _c: &Config,
) -> Vec<bool> {
    let l = {
        let l0 = answer.iter().take(training).filter(|&&e| e).count();
        let l1 = answer.iter().take(training).filter(|&&e| !e).count();
        l0.min(l1)
    };
    if l == 0 {
        (0..dataset.len()).map(|i| i % 2 == 0).collect()
    } else {
        let (d0, d1): (Vec<_>, Vec<_>) = dataset
            .iter()
            .map(|e| e.as_slice())
            .zip(answer.iter())
            .take(training)
            .partition(|(_, &b)| b);
        let mut d0: Vec<_> = d0.into_iter().map(|(x, _)| x).collect();
        let mut d1: Vec<_> = d1.into_iter().map(|(x, _)| x).collect();
        let mut _f = Factory::new();
        // let mut m0 = f.generate_from_ref(&d0[..l], k);
        // let mut m1 = f.generate_from_ref(&d1[..l], k);
        let mut is_updated = true;
        let mut pred: Vec<_> = vec![false; dataset.len() - training];
        while is_updated {
            is_updated = pred
                .iter_mut()
                .enumerate()
                .map(|(idx, p)| {
                    let chunk = &dataset[idx + training];
                    let f0 = d0.iter().map(|e| edlib_sys::global_dist(e, chunk)).min();
                    let f1 = d1.iter().map(|e| edlib_sys::global_dist(e, chunk)).min();
                    let next = f0 < f1;
                    let up = next != *p;
                    *p = next;
                    up
                })
                .fold(false, |p, q| p | q);
            d0 = dataset
                .iter()
                .zip(answer.iter().take(training).chain(pred.iter()))
                .filter_map(|(e, &b)| if b { Some(e.as_slice()) } else { None })
                .collect();
            d1 = dataset
                .iter()
                .zip(answer.iter().take(training).chain(pred.iter()))
                .filter_map(|(e, &b)| if !b { Some(e.as_slice()) } else { None })
                .collect();
            // m0 = f.generate_from_ref(&d0, k);
            // m1 = f.generate_from_ref(&d1, k);
        }
        answer.iter().take(training).copied().chain(pred).collect()
    }
}

