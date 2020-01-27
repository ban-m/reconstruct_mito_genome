extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
fn main() {
    let len = 50;
    let num_seq = 30;
    let test_num = 1000;
    let k = 6;
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12_218_993);
    let template1 = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let template2 = gen_sample::introduce_errors(&template1, &mut rng, 1, 1, 1);
    let template3 = gen_sample::introduce_errors(&template1, &mut rng, 1, 1, 1);
    assert_eq!(edlib_sys::global_dist(&template1, &template2), 3);
    assert_eq!(edlib_sys::global_dist(&template1, &template3), 3);
    let p = &gen_sample::PROFILE;
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, &mut rng, p))
        .collect();
    let data3: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template3, &mut rng, p))
        .collect();
    let weight = vec![1.; num_seq];
    let model2 = {
        let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
        DBGHMM::new_with_weight_prior(&data2, &weight, k)
    };
    let model3 = {
        let data3: Vec<_> = data3.iter().map(|e| e.as_slice()).collect();
        DBGHMM::new_with_weight_prior(&data3, &weight, k)
    };
    let tests: Vec<_> = (0..test_num)
        .map(|_| gen_sample::introduce_randomness(&template1, &mut rng, p))
        .collect();
    let tests2: Vec<_> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(0.5) {
                gen_sample::introduce_randomness(&template2, &mut rng, p)
            } else {
                gen_sample::introduce_randomness(&template3, &mut rng, p)
            }
        })
        .collect();
    println!("Entropy\tType\tSkew");
    eprintln!("{}\t{}", model2, model3);
    for (w2, w3) in tests.iter().map(|test| {
        let l2 = model2.forward(&test, &DEFAULT_CONFIG);
        let l3 = model3.forward(&test, &DEFAULT_CONFIG);
        eprintln!("TT:{:.3}\t{:.3}", l2, l3);
        as_weight(l2, l3)
    }) {
        let e = entropy(w2, w3);
        println!("{:.4}\tHMM\tNeutral", e);
    }
    for (w2, w3) in tests2.iter().map(|test| {
        let l2 = model2.forward(&test, &DEFAULT_CONFIG);
        let l3 = model3.forward(&test, &DEFAULT_CONFIG);
        as_weight(l2, l3)
    }) {
        let e = entropy(w2, w3);
        println!("{:.4}\tHMM\tSkew", e);
    }
    for (w1, w2) in tests.iter().map(|test| {
        let from2 = data2
            .iter()
            .map(|seq| edlib_sys::global_dist(seq, test))
            .min()
            .unwrap() as f64;
        let from3 = data3
            .iter()
            .map(|seq| edlib_sys::global_dist(seq, test))
            .min()
            .unwrap() as f64;
        let sum = from2.exp() + from3.exp();
        (from2.exp() / sum, from3.exp() / sum)
    }) {
        let e = entropy(w1, w2);
        println!("{:.4}\tAln\tNeutral", e);
    }
    for (w1, w2) in tests2.iter().map(|test| {
        let from2 = data2
            .iter()
            .map(|seq| edlib_sys::global_dist(seq, test))
            .min()
            .unwrap() as f64;
        let from3 = data3
            .iter()
            .map(|seq| edlib_sys::global_dist(seq, test))
            .min()
            .unwrap() as f64;
        let sum = from2.exp() + from3.exp();
        (from2.exp() / sum, from3.exp() / sum)
    }) {
        let e = entropy(w1, w2);
        println!("{:.4}\tAln\tSkew", e);
    }
}

fn entropy(x: f64, y: f64) -> f64 {
    if x < 0.00001 || y < 0.00001 {
        0.
    } else {
        -x * x.ln() - y * y.ln()
    }
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}
