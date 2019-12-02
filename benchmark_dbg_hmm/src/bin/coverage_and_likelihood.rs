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
    let mut rng: StdRng = SeedableRng::seed_from_u64(121_892);
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
    let k = 6;
    let mut f = Factory::new();
    debug!("{}", edlib_sys::global_dist(&template, &diff));
    println!("Coverage\tLikelihood\tType");
    for i in 1..num_seq {
        let m: Vec<_> = data[..i].iter().map(|e| e.as_slice()).collect();
        let w = vec![1.; i];
        let m = f.generate_with_weight(&m, &w, k);
        for _ in 0..100 {
            let q = introduce_randomness(&template, &mut rng, &PROFILE);
            let lk = m.forward(&q, &DEFAULT_CONFIG);
            println!("{}\t{:.4}\tTRUE", i, lk);
            let q = introduce_randomness(&diff, &mut rng, &PROFILE);
            let lk = m.forward(&q, &DEFAULT_CONFIG);
            println!("{}\t{:.4}\tFALSE", i, lk);
        }
    }
}
