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
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(i);
    let template = gen_sample::generate_seq(&mut rng, 150);
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    // for m in &model1 {
    //     eprintln!("{}", String::from_utf8_lossy(m));
    // }
    use std::time::Instant;
    let s = Instant::now();
    let model = POA::generate_vec(&model1);
    eprintln!("{:?}", Instant::now() - s);
    eprintln!("Model:{}\n{}", i, model);
    // eprintln!("Model:\n{:?}", model);
    // eprintln!("{}", String::from_utf8_lossy(&template));
    // eprintln!("{}", String::from_utf8_lossy(&model.consensus().unwrap()));
    // let (seq, model1) = model1.split_last().unwrap();
    // let model1 = model1.to_vec();
    // let model = POA::generate_vec(&model1);
    // eprintln!("{:?}", &model);
    // eprintln!("--------");
    // let model = model.add_with(seq, 1., &DEFAULT_CONFIG);
    // eprintln!("{:?}", model);
    // eprintln!("=========");
    // eprintln!("{:?}", model.remove_node());
    //}
}
