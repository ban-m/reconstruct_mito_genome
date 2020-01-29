extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let len = 150;
    let k = 6;
    let unit = 100;
    let rep: Vec<_> = (0..1_000).collect();
    let sum = rep
        .into_par_iter()
        .map(|i| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(1000 + i as u64);
            let template = gen_sample::generate_seq(&mut rng, len);
            let data1: Vec<_> = (0..unit)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
            let weight = vec![vec![1.; 10], vec![0.0001; unit - 10]].concat();
            let model1 = DBGHMM::new_with_weight(&data1, &weight, k);
            model1.node_num()
        })
        .sum::<usize>();
    eprintln!("{}", sum);
}
