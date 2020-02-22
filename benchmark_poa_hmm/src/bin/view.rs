extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::seq::SliceRandom;
use rand::{rngs::StdRng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    // let bases = b"ACTG";
    // let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    // let template: Vec<_> = (0..50)
    //     .filter_map(|_| bases.choose(&mut rng))
    //     .copied()
    //     .collect();
    // let model1: Vec<Vec<_>> = (0..13)
    //     .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
    //     .collect();
    let template = b"TACCTG".to_vec();
    let model1 = vec![
        template.clone(),
        b"TTACCTG".to_vec(),
        b"TAACCTG".to_vec(),
        b"TACCTG".to_vec(),
        // b"TACCTG".to_vec(),
    ];
    let model = POA::generate_vec(&model1);
    eprintln!("{}:TRUE", String::from_utf8_lossy(&template));
    eprintln!("{}:GEN", String::from_utf8_lossy(&model.gen().unwrap()));
    for seq in model1 {
        eprintln!("{}", String::from_utf8_lossy(&seq));
    }
    eprintln!("{}", model);
    eprintln!("{:?}", model);

    // let len = 20;
    // let mut rng: StdRng = SeedableRng::seed_from_u64(121332983);
    // let template: Vec<_> = generate_seq(&mut rng, len);
    // let p = Profile {
    //     sub: 0.01,
    //     ins: 0.01,
    //     del: 0.01,
    // };
    // println!("{}", String::from_utf8_lossy(&template));
    // let coverage = 60;
    // let data1: Vec<_> = (0..coverage)
    //     .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
    //     .collect();
    // let m1 = poa_hmm::POA::generate_vec(&data1);
    // println!("{:?}", m1);
}
