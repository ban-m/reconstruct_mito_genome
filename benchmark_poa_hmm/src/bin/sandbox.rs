extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::seq::SliceRandom;
use rand::{rngs::StdRng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

fn alignment<F>(xs: &[u8], ys: &[u8], del: f64, ins: f64, score: F) -> f64
where
    F: Fn(u8, u8) -> f64,
{
    let mut dp = vec![vec![0.; ys.len() + 1]; xs.len() + 1];
    for i in 0..xs.len() {
        dp[i + 1][0] = ins * (i + 1) as f64;
    }
    // for j in 0..ys.len() {
    //     dp[0][j + 1] = del * (j + 1) as f64;
    // }
    for i in 0..xs.len() {
        for j in 0..ys.len() {
            dp[i + 1][j + 1] = (dp[i][j] + score(xs[i], ys[j]))
                .max(dp[i + 1][j] + del)
                .max(dp[i][j + 1] + ins);
        }
    }
    *dp[xs.len()]
        .iter()
        .max_by(|a, b| a.partial_cmp(&b).unwrap())
        .unwrap()
}

fn main() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let mut model = POA::generate_vec(&model1);
    let s = std::time::Instant::now();
    eprintln!("{}", model);
    model.align(&model1[0], -1., -1., |x, y| if x == y { 1. } else { -1. });
    let e = std::time::Instant::now();
    eprintln!("{:?}", e - s);
    let q = gen_sample::generate_seq(&mut rng, 150);
    let q2 = gen_sample::generate_seq(&mut rng, 500);
    let s = std::time::Instant::now();
    let _ = alignment(&q, &q2, -1., -1., |x, y| if x == y { 1. } else { -1. });
    let e = std::time::Instant::now();
    eprintln!("{:?}", e - s);
}
