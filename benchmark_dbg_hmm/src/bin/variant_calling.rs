extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let k = 6;
    let mut rng: StdRng = SeedableRng::seed_from_u64(1213329);
    let chain_len = 40;
    let template: Vec<_> = (0..chain_len)
        .map(|_| generate_seq(&mut rng, len))
        .collect();
    let p = Profile {
        sub: 0.0003,
        ins: 0.0003,
        del: 0.0003,
    };
    let template2: Vec<_> = template
        .iter()
        .map(|t| introduce_randomness(t, &mut rng, &p))
        .collect();
    let dists: Vec<_> = template
        .iter()
        .zip(template2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .collect();
    let d = dists.iter().sum::<u32>();
    println!("{}", d);
    let coverage1 = 30;
    let coverage2 = 30;
    let mut gen = |ts: &Vec<Vec<u8>>| {
        ts.iter()
            .map(|t| introduce_randomness(&t, &mut rng, &PROFILE))
            .collect::<Vec<_>>()
    };
    let data1: Vec<_> = (0..coverage1).map(|_| gen(&template)).collect();
    let data2: Vec<_> = (0..coverage2).map(|_| gen(&template2)).collect();
    let data = {
        let mut data = vec![vec![]; template.len()];
        for d in data1.iter().chain(data2.iter()) {
            for (idx, c) in d.iter().enumerate() {
                data[idx].push(c.as_slice());
            }
        }
        data
    };
    let weight1 = vec![vec![1.2; coverage1], vec![0.2; coverage2]].concat();
    let weight2 = vec![vec![0.2; coverage1], vec![1.2; coverage2]].concat();
    let mut f = Factory::new();
    let mut buf = vec![];
    let m1: Vec<_> = data
        .iter()
        .map(|chunks| {
            let mut m = f.generate_with_weight(chunks, &weight1, k, &mut buf);
            m.check(&chunks, &DEFAULT_CONFIG, -120.);
            m
        })
        .collect();
    let m2: Vec<_> = data
        .iter()
        .map(|chunks| {
            let mut m = f.generate_with_weight(chunks, &weight2, k, &mut buf);
            m.check(&chunks, &DEFAULT_CONFIG, -120.);
            m
        })
        .collect();
    for (idx, (m1, m2)) in m1.iter().zip(m2.iter()).enumerate() {
        println!("{}\t{}\t{}", idx, m1, m2);
    }
    let test: Vec<_> = data1
        .iter()
        .chain(data2.iter())
        .cloned()
        .map(|chunks| last_decompose::ERead::new_with_lowseq(chunks, "0"))
        .collect();
    let models = vec![m1, m2];
    let test: Vec<_> = test
        .iter()
        .map(|read| {
            read.seq
                .iter()
                .map(|e| (e.unit(), e.bases().to_vec()))
                .collect()
        })
        .collect();
    let weight = last_decompose::variant_calling::variant_call(&models, &test, &DEFAULT_CONFIG);
    for (idx, (d, w)) in dists.iter().zip(weight.iter()).enumerate() {
        println!("{}\t{}\t{:.4}", idx, d, w);
    }
}
