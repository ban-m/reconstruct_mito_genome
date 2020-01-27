extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 50;
    let k = 6;
    let unit = 50;
    let _test = 10;
    let mut rng: StdRng = SeedableRng::seed_from_u64(1000);
    let template = gen_sample::generate_seq(&mut rng, len);
    let data1: Vec<_> = (0..unit * 2)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let template2 = template.clone();
    let data2: Vec<_> = (0..unit)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    // println!("{}", String::from_utf8_lossy(&template));
    // println!("{}", String::from_utf8_lossy(&template2));
    let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
    let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
    let weight = vec![1.; unit];
    let model1 = DBGHMM::new_with_weight_prior(&data1, &weight, k);
    let model2 = DBGHMM::new_with_weight_prior(&data2, &weight, k);
    println!("{}\n{:?}", model1, model1);
    println!("{}\n{:?}", model2, model2);
    // if model1.edge_num() != model2.edge_num() {
    //     let kmers = template
    //         .windows(k)
    //         .map(|kmer| kmer.to_vec())
    //         .collect::<HashSet<_>>()
    //         .len();
    //     println!("{}\n{}\n{}", kmers, model1, model2);
    // }
    // let tests: Vec<_> = (0..test)
    //     .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
    //     .collect();
    // let result: Vec<_> = tests
    //     .iter()
    //     .map(|q| {
    //         let l1 = model1.forward(q, &DEFAULT_CONFIG);
    //         let l2 = model2.forward(q, &DEFAULT_CONFIG);
    //         (l1, l2)
    //     })
    //     .collect();
    // for (l1, l2) in result {
    //     println!("{:.3}\t{:.3}", l1, l2);
    // }
}
