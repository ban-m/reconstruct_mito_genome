extern crate dbg_hmm;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let num_seq = 8;
    let k = 6;
    let mut rng: StdRng = SeedableRng::seed_from_u64(899_892);
    let template1: Vec<_> = generate_seq(&mut rng, len);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    for d in data1.iter() {
        eprintln!("{}", String::from_utf8_lossy(d));
    }
    let m1 = {
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        // DBGHMM::new_from_ref(&data1, k)
        let weight = vec![1.; num_seq];
        DBGHMM::new_with_weight(&data1, &weight, k)
    };
    // {
    //     let q = introduce_randomness(&template1, &mut rng, &PROFILE);
    //     let q = m1.forward(&q, &DEFAULT_CONFIG);
    let data = &data1[3];
    let r = m1.forward(&data, &DEFAULT_CONFIG);
    eprintln!("{}", r);
    let data = &data1[0];
    let r = m1.forward(&data, &DEFAULT_CONFIG);
    eprintln!("{}", r);
    // eprintln!("{:.4}\t{:.4}", q, r);
    //    }
}
