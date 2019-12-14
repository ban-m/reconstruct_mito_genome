extern crate dbg_hmm;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let num_seq = 50;
    let mut rng: StdRng = SeedableRng::seed_from_u64(899_892);
    let template1: Vec<_> = generate_seq(&mut rng, len);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let k = 6;
    let m1 = DBGHMM::new(&data1, k);
    for i in 0..10 {
        let q = introduce_randomness(&template1, &mut rng, &PROFILE);
        eprintln!("{}\t{:.4}", i % 2 + 1, m1.forward(&q, &DEFAULT_CONFIG),);
    }
}
