extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
#[macro_use]
extern crate log;
extern crate rand_xoshiro;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 50;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(121_892);
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let template: Vec<_> = generate_seq(&mut rng, len);
    let diff = introduce_randomness(&template, &mut rng, p);
    let data: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let data_diff: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&diff, &mut rng, &PROFILE))
        .collect();
    let k = 6;
    let mut f = Factory::new();
    debug!("{}", edlib_sys::global_dist(&template, &diff));
    println!("Coverage\tCorrect\tWrong");
    let tot = (1..num_seq)
        .map(|i| {
            let w = vec![1.; i];
            let m: Vec<_> = data[..i].iter().map(|e| e.as_slice()).collect();
            let m = f.generate_with_weight(&m, &w, k);
            let d: Vec<_> = data_diff[..i].iter().map(|e| e.as_slice()).collect();
            let d = f.generate_with_weight(&d, &w, k);
            let correct = (0..100)
                .map(|e| {
                    let q = introduce_randomness(&template, &mut rng, &PROFILE);
                    let lk = m.forward(&q, &DEFAULT_CONFIG);
                    let lkd = d.forward(&q, &DEFAULT_CONFIG);
                    println!("{}\t{}\t{}", i, lk, lkd);
                    lk > lkd
                })
                .filter(|&e| e)
                .count() as f64
                / 100.;
            debug!("{}\t{}", i, correct);
            correct
        })
        .sum::<f64>()
        / (num_seq - 1) as f64;
    debug!("Overall:{}", tot);
}
