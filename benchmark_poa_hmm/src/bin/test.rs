extern crate edlib_sys;
extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::gen_sample::*;
use poa_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let len = 150;
    let test = 1000;
    for i in 0..100 {
        let mut rng: StdRng = SeedableRng::seed_from_u64(i as u64);
        let coverage = 30;
        let template: Vec<_> = generate_seq(&mut rng, len);
        let data1: Vec<_> = (0..coverage)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let score = |x, y| if x == y { 3 } else { -4 };
        let param = (-6, -6, &score);
        let m1 = poa_hmm::POA::default().update(&data1, &vec![0.2; coverage], param);
        //eprintln!("{}", m1);
        let m2 = poa_hmm::POA::default().update(&data1, &vec![0.4; coverage], param);
        // eprintln!("{}", m2);
        // eprintln!("-----");
        let tests: Vec<_> = (0..test)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let count = tests
            .par_iter()
            .filter(|t| m1.forward(t, &DEFAULT_CONFIG) >= m2.forward(t, &DEFAULT_CONFIG))
            .count();
        eprintln!("{}/{}", count, test);
    }
}
