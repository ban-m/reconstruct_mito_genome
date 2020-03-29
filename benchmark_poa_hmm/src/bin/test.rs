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
    let test = 1000;
    let time = 4;
    let coverage = 5;
    for i in 0..100 {
        let mut rng: StdRng = SeedableRng::seed_from_u64(i as u64);
        let template: Vec<_> = generate_seq(&mut rng, len);
        // let template_sub = introduce_errors(&template, &mut rng, 1, 0, 0);
        // let template_del = introduce_errors(&template, &mut rng, 0, 1, 0);
        // let template_ins = introduce_errors(&template, &mut rng, 0, 0, 1);
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
        // let data_sub: Vec<_> = (0..coverage)
        //     .map(|_| introduce_randomness(&template_sub, &mut rng, &PROFILE))
        //     .collect();
        // let data_del: Vec<_> = (0..coverage)
        //     .map(|_| introduce_randomness(&template_del, &mut rng, &PROFILE))
        //     .collect();
        // let data_ins: Vec<_> = (0..coverage)
        //     .map(|_| introduce_randomness(&template_ins, &mut rng, &PROFILE))
        //     .collect();
        let score = |x, y| if x == y { 3 } else { -4 };
        let param = (-6, -6, &score);
        let prior = poa_hmm::POA::generate(&dataset, &vec![1.; coverage * (time + 1)], param);
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
        let m = prior.clone().update(&data1, &vec![1.; coverage], param, 1);
        let m_high = prior
            .clone()
            .update(&data2, &vec![1.; coverage * time], param, 1);
        // let tests: Vec<_> = (0..test)
        //     .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        //     .collect();
        let result: Vec<_> = data1
            .par_iter()
            .map(|q| {
                let weight_ratio = 20. / m.weight();
                let prior = (1. + weight_ratio).ln();
                let m1 = m.forward(q, &DEFAULT_CONFIG) - prior;
                let weight_ratio = 20. / m_high.weight();
                let prior = (1. + weight_ratio).ln();
                let m2 = m_high.forward(q, &DEFAULT_CONFIG) - prior;
                (m1, m2, m1 - m2)
            })
            .collect();
        for (m1, m2, diff) in result {
            println!("{:.3}\t{:.3}\t{:.3}", m1, m2, diff);
        }
        // let m_sub = poa_hmm::POA::from_vec(&data_sub, &vec![1.; coverage], param);
        // let m_del = poa_hmm::POA::from_vec(&data_del, &vec![1.; coverage], param);
        // let m_ins = poa_hmm::POA::from_vec(&data_ins, &vec![1.; coverage], param);
        // eprintln!(" :{}\nS:{}\nD:{}\nI:{}", m, m_sub, m_del, m_ins);
        // let diffs: Vec<_> = tests
        //     .par_iter()
        //     .map(|q| {
        //         let normal = m.forward(q, &DEFAULT_CONFIG);
        //         let sub = m_sub.forward(q, &DEFAULT_CONFIG);
        //         let del = m_del.forward(q, &DEFAULT_CONFIG);
        //         let ins = m_ins.forward(q, &DEFAULT_CONFIG);
        //         (normal - sub, normal - del, normal - ins)
        //     })
        //     .collect();
        // for (sub, del, ins) in diffs {
        //     eprintln!("DIFF\t{:.2}\t{:.2}\t{:.2}", sub, del, ins);
        // }
    }
}
