extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::seq::SliceRandom;
use rand::{rngs::StdRng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
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
    let model = POA::generate_vec(&model1);
    println!("{}/{:?}", model, model);
}
