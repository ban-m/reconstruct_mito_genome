extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
#[macro_use]
extern crate log;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 50;
    let mut rng: StdRng = SeedableRng::seed_from_u64(12_121_899_892);
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
    let full: Vec<_> = data
        .iter()
        .chain(data_diff.iter())
        .map(|e| e.as_slice())
        .collect();
    let ws: Vec<_> = vec![1.; num_seq * 2];
    let full = f.generate_with_weight(&full, &ws, k);
    debug!("{}", edlib_sys::global_dist(&template, &diff));
    println!("Coverage\tCorrect\tWrong\tFull");
    for i in 1..num_seq {
        let w = vec![1.; i];
        let m: Vec<_> = data[..i].iter().map(|e| e.as_slice()).collect();
        let m = f.generate_with_weight(&m, &w, k);
        let d: Vec<_> = data_diff[..i].iter().map(|e| e.as_slice()).collect();
        let d = f.generate_with_weight(&d, &w, k);
        for _ in 0..100 {
            let q = introduce_randomness(&template, &mut rng, &PROFILE);
            let lk = m.forward(&q, &DEFAULT_CONFIG);
            let lkd = d.forward(&q, &DEFAULT_CONFIG);
            let lkf = full.forward(&q, &DEFAULT_CONFIG);
            println!("{}\t{}\t{}\t{}", i, lk, lkd, lkf);
        }
    }
}
