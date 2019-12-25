extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::gen_sample;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
// use std::time::Instant;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let len = 150;
    let k = 6;
    let num_seq = 100;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let c = dbg_hmm::STRICT_CONFIG;
    // let c = dbg_hmm::DEFAULT_CONFIG;
    let s = &gen_sample::PROFILE;
    for i in 0..10 {
        // let start = Instant::now();
        let (acc, d) = benchmark(p, s, i as u64 * 849, k, len, num_seq, &c);
        debug!("Acc:{:.4}\t{}", acc, d);
        // let ela = Instant::now() - start;
        // let millis = ela.as_millis() as f64 / num_seq as f64;
        // info!("End {}, {:?}({:.4}/sample)", num_seq, ela, millis);
    }
}
use dbg_hmm::Config;
fn benchmark(
    p: &gen_sample::Profile,
    s: &gen_sample::Profile,
    seed: u64,
    k: usize,
    len: usize,
    num_seq: usize,
    config: &Config,
) -> (f64, u32) {
    let cluster_num = 2;
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(138222324 + seed);
    let template = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let templates: Vec<_> = (0..cluster_num)
        .map(|_| dbg_hmm::gen_sample::introduce_randomness(&template, &mut rng, &p))
        .collect();
    let answers: Vec<_> = (0..num_seq)
        .map(|_| rng.gen_range(0, cluster_num))
        .collect();
    let dataset: Vec<Vec<_>> = answers
        .iter()
        .map(|&e| gen_sample::introduce_randomness(&templates[e], &mut rng, s))
        .collect();
    let mut d = 0;
    for (i, seqi) in templates.iter().enumerate() {
        for (j, seqj) in templates.iter().enumerate().skip(i + 1) {
            let dist = edlib_sys::global_dist(seqi, seqj);
            debug!("{}-{}:{}", i, j, dist);
            d += dist;
        }
    }
    let data: Vec<_> = dataset.iter().map(|e| e.as_slice()).collect();
    let pred = last_decompose::unit_clustering(&data, k, cluster_num + 2, config);
    let mut result = vec![vec![0; cluster_num + 2]; cluster_num + 2];
    for (&pred, &ans) in pred.iter().zip(answers.iter()) {
        result[pred as usize][ans as usize] += 1;
    }
    for (idx, preds) in result.iter().enumerate() {
        let tp = preds[idx];
        let tot = preds.iter().sum::<u32>();
        print!("Predicted as {}:", idx);
        for ans in preds {
            print!("{}\t", ans);
        }
        println!("Total:{:.4}", tp as f64 / tot as f64);
    }
    let acc = pred
        .iter()
        .zip(answers)
        .filter(|&(&p, a)| p == a as u8)
        .count() as f64
        / num_seq as f64;
    let cluster_num = cluster_num as u32;
    (acc, d / (cluster_num + 1) / cluster_num)
}
