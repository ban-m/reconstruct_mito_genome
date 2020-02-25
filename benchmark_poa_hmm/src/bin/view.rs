extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    // let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(121212);
    // let template: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
    // let model = vec![
    //     b"CAGTGTCAGTGCTAGCT".to_vec(),
    //     b"CGTGTCAGTGCTAGCT".to_vec(),
    //     b"CAGGTCAGTGCTAGCT".to_vec(),
    // ];
    // eprintln!(" :{}", String::from_utf8_lossy(&template));
    // for m in &model {
    //     eprintln!("1:{}", String::from_utf8_lossy(m));
    // }
    // let model = POA::generate_vec(&model);
    // eprintln!("Model1:\n{:?}", model);
    // eprintln!(" :{}", String::from_utf8_lossy(&template));
    // eprintln!(" :{}", String::from_utf8_lossy(&model.consensus().unwrap()));

    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212321);
    let template = gen_sample::generate_seq(&mut rng, 150);
    let model1: Vec<Vec<_>> = (0..25)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    for m in &model1 {
        eprintln!("1:{}", String::from_utf8_lossy(m));
    }
    let model1 = POA::generate_vec(&model1);
    eprintln!("Model:\n{:?}", model1);
    eprintln!(" :{}", String::from_utf8_lossy(&template));
    eprintln!(
        " :{}",
        String::from_utf8_lossy(&model1.consensus().unwrap())
    );
}
