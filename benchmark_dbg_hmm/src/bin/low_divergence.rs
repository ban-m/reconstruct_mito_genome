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
    let unit = 25;
    let rep: Vec<_> = (0..20000).collect();
    let frac = 40;
    let result: Vec<_> = rep
        .into_par_iter()
        .map(|i| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(121332983 + i);
            let template = generate_seq(&mut rng, len);
            let data1: Vec<_> = (0..unit)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let template2 = match i % frac {
                0 => introduce_errors(&template, &mut rng, 1, 0, 0),
                1 => introduce_errors(&template, &mut rng, 0, 1, 0),
                2 => introduce_errors(&template, &mut rng, 0, 0, 1),
                _ => template.clone(),
            };
            let data2: Vec<_> = (0..unit)
                .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
                .collect();
            let data: Vec<_> = data1
                .iter()
                .chain(data2.iter())
                .map(|e| e.as_slice())
                .collect();
            let mut weight1 = vec![];
            let mut weight2 = vec![];
            for _ in 0..unit * 2 {
                use rand::Rng;
                let w = rng.gen_range(0.4, 0.6);
                weight1.push(w);
                weight2.push(1. - w);
            }
            let m1 = DBGHMM::new_with_weight(&data, &weight1, k);
            let m2 = DBGHMM::new_with_weight(&data, &weight2, k);
            let test: Vec<_> = (0..20)
                .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                .collect();
            let lkdiff: Vec<_> = test
                .par_iter()
                .map(|d| {
                    let f1 = m1.forward(&d, &dbg_hmm::DEFAULT_CONFIG);
                    let f2 = m2.forward(&d, &dbg_hmm::DEFAULT_CONFIG);
                    (f1 - f2).abs()
                })
                .collect();
            let len = lkdiff.len() as f64;
            let mean = lkdiff.iter().sum::<f64>() / len;
            let sd = lkdiff.iter().map(|l| (l - mean).powi(2i32)).sum::<f64>() / len;
            let t = match i % frac {
                0 => 's',
                1 => 'd',
                2 => 'i',
                _ => 'n',
            };
            format!("{}\t{}\t{}", t, mean, sd)
        })
        .collect();
    println!("Type\tMean\tSD");
    for line in result {
        println!("{}", line);
    }
}
