#![feature(test)]
extern crate packed_simd;
extern crate poa_hmm;
extern crate rand;
extern crate rand_xoshiro;
extern crate test;
use poa_hmm::gen_sample::*;
use poa_hmm::*;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use test::bench::Bencher;

fn alignment_float<F>(xs: &[u8], ys: &[u8], del: f64, ins: f64, score: F) -> f64
where
    F: Fn(u8, u8) -> f64,
{
    let mut dp = vec![vec![0.; ys.len() + 1]; xs.len() + 1];
    for i in 0..xs.len() {
        dp[i + 1][0] = ins * (i + 1) as f64;
    }
    for j in 0..ys.len() {
        dp[0][j + 1] = del * (j + 1) as f64;
    }
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

fn alignment_int<F>(xs: &[u8], ys: &[u8], del: i32, ins: i32, score: F) -> i32
where
    F: Fn(u8, u8) -> i32,
{
    let (row, column) = (xs.len() + 1, ys.len() + 1);
    let mut dp = vec![0; row * column];
    for i in 0..xs.len() {
        dp[(i + 1) * column] = ins * (i + 1) as i32;
    }
    for j in 0..ys.len() {
        dp[j + 1] = del * (j + 1) as i32;
    }
    for i in 0..xs.len() {
        let current_row = (i + 1) * column;
        let previous_row = current_row - column;
        for j in 0..ys.len() {
            let current_pos = current_row + j;
            let previous_pos = previous_row + j;
            dp[current_pos + 1] = (dp[previous_pos] + score(xs[i], ys[j]))
                .max(dp[current_pos] + del)
                .max(dp[previous_pos + 1] + ins);
        }
    }
    *dp[xs.len() * column..].iter().max().unwrap()
}

#[bench]
fn alignment_i32(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template1: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let template2: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    //let score = |x, y| SCORE[(x == y) as usize];
    let score = |x, y| if x == y { 3 } else { -1 };
    b.iter(|| alignment_int(&template1, &template2, -2, -2, score))
}

#[bench]
fn alignment_f64(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template1: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let template2: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let score = |x, y| if x == y { 3. } else { -1. };
    b.iter(|| alignment_float(&template1, &template2, -2., -2., score))
}

#[bench]
fn align_150(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let model = POA::generate_vec(&model1);
    b.iter(|| {
        model
            .clone()
            .align(&model1[0], (-1, -1, &|x, y| if x == y { 1 } else { -1 }))
    });
}

#[bench]
fn add_150(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let model = POA::generate_vec(&model1);
    b.iter(|| {
        model.clone().add(
            &model1[0],
            1.,
            (-1, -1, &|x, y| if x == y { 1 } else { -1 }),
        )
    });
}

#[bench]
fn create_150(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    b.iter(|| POA::generate_vec(&model1));
}

#[bench]
fn forward_150(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let query = introduce_randomness(&template, &mut rng, &PROFILE);
    let model = POA::generate_vec(&model1);
    eprintln!("{}", model);
    b.iter(|| model.forward(&query, &DEFAULT_CONFIG));
}
