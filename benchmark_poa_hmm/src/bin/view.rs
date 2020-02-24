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

    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212);
    let template = b"AGGACCGTCTCTTCTCCAGTCTACG";
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(template, &mut rng, &PROFILE))
        .collect();
    // let model2: Vec<Vec<_>> = (0..10)
    //     .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
    //     .collect();
    for m in &model1 {
        eprintln!("1:{}", String::from_utf8_lossy(m));
    }
    // eprintln!(" :{}", String::from_utf8_lossy(&template1));
    // for m in &model2 {
    //     eprintln!("2:{}", String::from_utf8_lossy(m));
    // }
    // let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
    // let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
    // eprintln!("1:{}", String::from_utf8_lossy(&test1));
    // eprintln!("2:{}", String::from_utf8_lossy(&test2));
    let model1 = POA::generate_vec(&model1);
    // eprintln!("Model1:\n{:?}", model1);
    // let model2 = POA::generate_vec(&model2);
    eprintln!("Model:\n{:?}", model1);
    eprintln!(" :{}", String::from_utf8_lossy(template));
    eprintln!(
        " :{}",
        String::from_utf8_lossy(&model1.consensus().unwrap())
    );
}
