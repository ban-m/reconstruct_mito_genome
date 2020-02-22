extern crate edlib_sys;
extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::gen_sample::*;
use poa_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let len = 150;
    let mut rng: StdRng = SeedableRng::seed_from_u64(12133);
    let template: Vec<_> = generate_seq(&mut rng, len);
    let p = Profile {
        sub: 0.003,
        ins: 0.003,
        del: 0.003,
    };
    let template2 = introduce_randomness(&template, &mut rng, &p);
    let d = edlib_sys::global_dist(&template, &template2);
    println!("Dist:{}", d);
    let coverage = 20;
    let data1: Vec<_> = (0..coverage)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let data2: Vec<_> = (0..coverage)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    let m1 = poa_hmm::POA::generate_vec(&data1);
    let m2 = poa_hmm::POA::generate_vec(&data2);
    println!("{}\n--------------\n{}", m1, m2);
    let num = 200;
    let test: Vec<_> = (0..2 * num)
        .map(|i| {
            if i % 2 == 0 {
                introduce_randomness(&template, &mut rng, &PROFILE)
            } else {
                introduce_randomness(&template2, &mut rng, &PROFILE)
            }
        })
        .collect();
    let correct1 = test
        .par_iter()
        .enumerate()
        .filter(|(idx, _)| idx % 2 == 0)
        .filter(|&(_, q)| {
            let f1 = m1.forward(q, &DEFAULT_CONFIG);
            let f2 = m2.forward(q, &DEFAULT_CONFIG);
            eprintln!("{}\t{}\t{}", 0, f1, f2);
            f1 > f2
        })
        .count();
    let correct2 = test
        .par_iter()
        .enumerate()
        .filter(|(idx, _)| idx % 2 == 1)
        .filter(|&(_, q)| {
            let f1 = m1.forward(q, &DEFAULT_CONFIG);
            let f2 = m2.forward(q, &DEFAULT_CONFIG);
            eprintln!("{}\t{}\t{}", 1, f1, f2);
            f1 < f2
        })
        .count();
    println!("Fin:{}/{}", correct1, correct2);
}
