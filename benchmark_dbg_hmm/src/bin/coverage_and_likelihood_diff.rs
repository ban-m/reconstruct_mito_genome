extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
#[macro_use]
extern crate log;
extern crate rand_xoshiro;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 200;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(121_899_200);
    let p = &gen_sample::Profile {
        sub: 0.004,
        ins: 0.004,
        del: 0.004,
    };
    let template: Vec<_> = generate_seq(&mut rng, len);
    let diff = introduce_randomness(&template, &mut rng, p);
    //let diff = introduce_errors(&template[..4], &mut rng, 0, 1, 0);
    //let diff: Vec<_> = diff.iter().chain(template[4..].iter()).copied().collect();
    let data: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let data_diff: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&diff, &mut rng, &PROFILE))
        .collect();
    let k = 6;
    let mut f = Factory::new();
    debug!("Dist:{}", edlib_sys::global_dist(&template, &diff));
    println!("Coverage\tCorrect\tWrong");
    let tot = (1..num_seq)
        .map(|i| {
            let w = vec![1.; i];
            let m: Vec<_> = data[..i].iter().map(|e| e.as_slice()).collect();
            let m = f.generate_with_weight_prior(&m, &w, k, &mut vec![]);
            let d: Vec<_> = data_diff[..i].iter().map(|e| e.as_slice()).collect();
            let d = f.generate_with_weight_prior(&d, &w, k, &mut vec![]);
            let correct = (0..100)
                .map(|_| {
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
