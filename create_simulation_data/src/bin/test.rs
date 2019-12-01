extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use create_simulation_data::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
// use rayon::prelude::*;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let test_num = 300;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    use std::time::Instant;
    let seed = 11920981;
    let coverage = 6;
    let s = Instant::now();
    let (hmm, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
    eprintln!("{:?}", Instant::now() - s);
    info!("Cov:{}", coverage);
    info!("Acc:{:.4}", hmm);
    info!("Dist:{} out of {}", dist, chain_len);
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
    let seed = 10034374 + seed;
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
    template1
        .iter()
        .zip(template2.iter())
        .enumerate()
        .for_each(|(idx, (t1, t2))| debug!("{}\t{}", idx, edlib_sys::global_dist(t1, t2)));
    let prob_0 = 0.5;
    let (dataset, label, answer, border) =
        generate_dataset(&template1, &template2, coverage, test_num, &mut rng, prob_0);
    // let pred = align_solve(&dataset, &label, border);
    // use rand::Rng;
    // let correct = pred
    //     .iter()
    //     .zip(answer.iter())
    //     .filter(|(p, a)| a == p)
    //     .count();
    // debug!("Initial:{:.4}", correct as f64 / test_num as f64);
    // {
    //     let m0: Vec<_> = label
    //         .iter()
    //         .chain(answer.iter())
    //         .map(|&e| if e == 0 { 1. } else { 0. })
    //         .collect();
    //     let w0 = m0.iter().sum::<f64>() / m0.len() as f64;
    //     let m1: Vec<_> = label
    //         .iter()
    //         .chain(answer.iter())
    //         .map(|&e| if e == 1 { 1. } else { 0. })
    //         .collect();
    //     let m0 = construct_with_weights(&dataset, &m0, k);
    //     let m1 = construct_with_weights(&dataset, &m1, k);
    //     let obj = calc_lk(&dataset[coverage..], &m0, &m1, w0, 1. - w0);
    //     debug!("OBJLK:{:.4}", obj);
    // }
    let em_pred = em_solve(&dataset, &label, border, k, &answer);
    let pos = answer.iter().filter(|&&e| e == 0).count();
    let neg = answer.len() - pos;
    let tp = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|&(&p, &a)| a == p && a == 0)
        .count();
    let tn = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|&(&p, &a)| a == p && a == 1)
        .count();
    debug!("Pos\tNeg\tTP\tTN");
    debug!("{}\t{}\t{}\t{}", pos, neg, tp, tn);
    let correct = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|(p, a)| a == p)
        .count();
    let hmm = correct as f64 / test_num as f64;
    (hmm, dist)
}
