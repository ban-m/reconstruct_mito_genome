extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::gen_sample;
use last_decompose::{clustering, likelihood_of_assignments, ERead};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (test_num, coverage, prob) = if args.len() > 3 {
        let tn = args[1].parse::<usize>().unwrap();
        let cov = args[2].parse::<usize>().unwrap();
        let prob = args[3].parse::<f64>().unwrap();
        (tn, cov, prob)
    } else {
        (200, 0, 0.5)
    };
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    use std::time::Instant;
    let seed = 11920981;
    let s = Instant::now();
    let (hmm, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len, prob);
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
    prob: f64,
) -> (f64, u32) {
    let seed = 1003437 + seed;
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
    let (dataset, label, answer, _border) = create_simulation_data::generate_dataset(
        &template1, &template2, coverage, test_num, &mut rng, prob,
    );
    let contigs = vec![chain_len];
    let data: Vec<_> = dataset
        .into_iter()
        .enumerate()
        .map(|(idx, e)| {
            let id = format!("{}", idx);
            ERead::new_with_lowseq(e, &id)
        })
        .collect();
    {
        let answer: Vec<_> = label.iter().chain(answer.iter()).copied().collect();
        let objlk = likelihood_of_assignments(&data, &answer, k, 2, &contigs);
        debug!("ObjLK:{}", objlk);
    }
    let forbidden = vec![vec![]; data.len()];
    let em_pred = clustering(&data, &label, &forbidden, k, 2, &contigs);
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
