extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    // for i in 0..100 {
    //     let i = i as u64;
    let i = 93;
    let coverage = 2;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(i);
    let template = gen_sample::generate_seq(&mut rng, 150);
    let model1: Vec<Vec<_>> = (0..coverage)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let model1: Vec<_> = model1.iter().map(|e| e.as_slice()).collect();
    use std::time::Instant;
    let score = |x, y| if x == y { 3 } else { -4 };
    let ws = vec![1.; coverage];
    // let s = Instant::now();
    // let model = POA::generate_w_param(&model1, &ws, -6, -6, &score);
    // eprintln!("{:?}", Instant::now() - s);
    // eprintln!("Model:{}\n{}", i, model);
    // let s = Instant::now();
    let mut model = POA::generate_w_param_simd(&model1, &ws, -6, -6, &score);
    // eprintln!("{:?}", Instant::now() - s);
    // eprintln!("Model:{}\n{}", i, model);
    let test = introduce_randomness(&template, &mut rng, &PROFILE);
    let (score, ops) = model.align_simd(&test, -6, -6, &score);
    let (q, g) = model.view(&test, &ops);
    eprintln!("{}", score);
    eprintln!("{}", q);
    eprintln!("{}", g);
}
