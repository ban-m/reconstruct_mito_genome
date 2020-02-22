#![feature(test)]
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
        let mut m = model.clone();
        m.align(&model1[0], -1., -1., |x, y| if x == y { 1. } else { -1. })
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
    b.iter(|| model.clone().add(&model1[0], 1.));
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
