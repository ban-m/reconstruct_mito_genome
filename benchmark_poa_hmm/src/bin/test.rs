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
    let mut rng: StdRng = SeedableRng::seed_from_u64(12133);
    let coverage = 30;
    for _ in 0..20 {
        let template: Vec<_> = generate_seq(&mut rng, len);
        let data1: Vec<_> = (0..coverage)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let score = |x, y| if x == y { 3 } else { -5 };
        let m1 = poa_hmm::POA::default();
        let m2 = poa_hmm::POA::default();
        let m1 = m1.update_s(&data1, &vec![0.4; coverage], (-4, -4, &score), 0);
        let m2 = m2.update_s(&data1, &vec![0.2; coverage], (-4, -4, &score), 0);
        eprintln!("{}\t{}", m1, m2);
        let num = 1000;
        let test: Vec<_> = (0..num)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let lks: Vec<_> = test
            .par_iter()
            .map(|q| {
                let f1 = m1.forward(q, &DEFAULT_CONFIG);
                let f2 = m2.forward(q, &DEFAULT_CONFIG);
                (f1, f2)
            })
            .collect();
        let count = lks
            .iter()
            .filter(|&(lk1, lk2)| {
                //println!("{}\t{}", lk1, lk2);
                lk1 < lk2
            })
            .count();
        eprintln!("{}/{}", count, num);
    }
}
