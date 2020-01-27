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
    let unit = 100;
    let test = 1000;
    for i in 0..100 {
        let mut rng: StdRng = SeedableRng::seed_from_u64(i);
        let template = generate_seq(&mut rng, len);
        let data1: Vec<_> = (0..unit)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let data2: Vec<_> = (0..unit)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
        let weight = vec![1.; unit];
        let model1 = DBGHMM::new_with_weight_prior(&data1, &weight, k);
        let model2 = DBGHMM::new_with_weight_prior(&data2, &weight, k);
        // if model1.edge_num() != model2.edge_num() {
        //     let kmers = template
        //         .windows(k)
        //         .map(|kmer| kmer.to_vec())
        //         .collect::<HashSet<_>>()
        //         .len();
        //     eprintln!("{}\t{}\n{}\n{}", i, kmers, model1, model2);
        // }
        // COMMM
        // eprintln!("{}", String::from_utf8_lossy(&template));
        // model1.dump(&template);
        // eprintln!("======");
        // model2.dump(&template);
        let tests: Vec<_> = (0..test)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let result: Vec<_> = tests
            .par_iter()
            .map(|q| {
                let l1 = model1.forward(q, &DEFAULT_CONFIG);
                let l2 = model2.forward(q, &DEFAULT_CONFIG);
                l1 - l2
            })
            .collect();
        for r in result {
            println!("{}\t{}", i, r);
        }
    }
}
