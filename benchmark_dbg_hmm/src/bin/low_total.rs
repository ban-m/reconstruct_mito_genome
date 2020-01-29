extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let len = 50;
    let k = 6;
    let unit = 125;
    let rep: Vec<_> = (0..20000).collect();
    let test_num = 100;
    let result: Vec<_> = rep
        .into_par_iter()
        .map(|i| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(121332983 + i);
            let template = generate_seq(&mut rng, len);
            use std::collections::HashSet;
            let data1: Vec<_> = (0..unit)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let data2: Vec<_> = (0..unit)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let kmers = template
                .windows(k)
                .map(|kmer| kmer.to_vec())
                .collect::<HashSet<_>>()
                .len();
            let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
            let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
            let weight = vec![1.; unit];
            let model1 = DBGHMM::new_with_weight(&data1, &weight, k);
            let model2 = DBGHMM::new_with_weight(&data2, &weight, k);
            eprintln!("{}\n{}\n{}", kmers, model1, model2);
            eprintln!("=====");
            let test: Vec<_> = (0..test_num)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let lks: Vec<_> = test
                .par_iter()
                .map(|d| {
                    let l1 = model1.forward(&d, &dbg_hmm::DEFAULT_CONFIG);
                    let l2 = model2.forward(&d, &dbg_hmm::DEFAULT_CONFIG);
                    l1 - l2
                })
                .collect();
            let max_f = |x, &y| if x < y { y } else { x };
            let min_f = |x, &y| if x < y { x } else { y };
            let max = lks.iter().fold(-10000., max_f);
            let min = lks.iter().fold(10000., min_f);
            let mean = lks.iter().sum::<f64>() / lks.len() as f64;
            let sumsq = lks.iter().map(|&e| e * e).sum::<f64>() / lks.len() as f64;
            let var = sumsq - mean * mean;
            format!("{}\t{}\t{}\t{}", max, min, mean, var)
        })
        .collect();
    println!("Max\tMin\tMean\tVar");
    for line in result {
        println!("{}", line);
    }
}
