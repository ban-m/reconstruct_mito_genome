extern crate edlib_sys;
extern crate poa_hmm;
extern crate rand;
use poa_hmm::gen_sample::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let time = 4;
    let coverage = 5;
    for i in 0..100 {
        let mut rng: StdRng = SeedableRng::seed_from_u64(i as u64);
        let template: Vec<_> = generate_seq(&mut rng, len);
        let data1: Vec<_> = (0..coverage)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let data2: Vec<_> = (0..coverage * time)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let dataset: Vec<_> = data1
            .iter()
            .chain(data2.iter())
            .map(|e| e.as_slice())
            .collect();
        let score = |x, y| if x == y { 3 } else { -4 };
        let param = (-6, -6, &score);
        let prior = poa_hmm::POA::generate(&dataset, &vec![1.; coverage * (time + 1)], param);
        eprintln!("{}", prior);
    }
}
