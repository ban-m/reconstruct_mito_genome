#![feature(test)]
extern crate dbg_hmm;
extern crate test;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng}; // 0.2%, 0.65%, 0.65%.
const SUB: f64 = 0.002;
const DEL: f64 = 0.0065;
const IN: f64 = 0.0065;
use dbg_hmm::DEFAULT_CONFIG;
use dbg_hmm::*;
use test::Bencher;
#[bench]
fn determine(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: StdRng = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..30)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..20)
        .map(|_| introduce_randomness(&template, &mut rng))
        .collect();
    let k = 7;
    let model1 = DBGHMM::new(&model1, k);
    b.iter(|| test::black_box(model1.forward(&template, &DEFAULT_CONFIG)));
}

#[bench]
fn new(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: StdRng = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..30)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..20)
        .map(|_| introduce_randomness(&template, &mut rng))
        .collect();
    let k = 7;
    b.iter(|| test::black_box(DBGHMM::new(&model1, k)));
}

#[bench]
fn new2(b: &mut Bencher) {
    let bases = b"ACTG";
    let mut rng: StdRng = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..30)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..20)
        .map(|_| introduce_randomness(&template, &mut rng))
        .collect();
    let k = 7;
    let mut f = Factory::new();
    b.iter(|| test::black_box(f.generate(&model1, k)));
}

enum Op {
    Match,
    MisMatch,
    Del,
    In,
}
impl Op {
    fn weight(&self) -> f64 {
        match self {
            Op::Match => 1. - SUB - DEL - IN,
            Op::MisMatch => SUB,
            Op::Del => DEL,
            Op::In => IN,
        }
    }
}
const OPERATIONS: [Op; 4] = [Op::Match, Op::MisMatch, Op::Del, Op::In];
fn introduce_randomness<T: rand::Rng>(seq: &[u8], rng: &mut T) -> Vec<u8> {
    let mut res = vec![];
    let mut remainings: Vec<_> = seq.iter().copied().rev().collect();
    while !remainings.is_empty() {
        match OPERATIONS.choose_weighted(rng, Op::weight).unwrap() {
            &Op::Match => res.push(remainings.pop().unwrap()),
            &Op::MisMatch => res.push(choose_base(rng, remainings.pop().unwrap())),
            &Op::In => res.push(random_base(rng)),
            &Op::Del => {
                remainings.pop().unwrap();
            }
        }
    }
    res
}
fn choose_base<T: rand::Rng>(rng: &mut T, base: u8) -> u8 {
    let bases: Vec<u8> = b"ATCG".iter().filter(|&&e| e != base).copied().collect();
    *bases.choose_weighted(rng, |_| 1. / 3.).unwrap()
}
fn random_base<T: rand::Rng>(rng: &mut T) -> u8 {
    *b"ATGC".choose_weighted(rng, |_| 1. / 4.).unwrap()
}
