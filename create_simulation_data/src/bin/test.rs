extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use create_simulation_data::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1000);
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let chain_len = 40;
    let len = 150;
    let coverage = 0;
    let test_num = 30;
    let _class = 3;
    let template1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let template2: Vec<_> = template1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let (dataset, _, answer, _) = generate_mul_data(
        &[template1.clone(), template2.clone()],
        coverage,
        test_num,
        &mut rng,
        &[0.5, 0.5],
    );
    let data: Vec<Vec<_>> = dataset
        .into_iter()
        .map(|read| read.into_iter().enumerate().collect())
        .collect();
    let models: Vec<Vec<poa_hmm::POA>> = (0..2)
        .map(|cl| {
            let cl = cl as u8;
            let mut chunks = vec![vec![]; chain_len];
            for (seq, _) in data.iter().zip(answer.iter()).filter(|&(_, &b)| b == cl) {
                for &(idx, ref chunk) in seq.iter() {
                    chunks[idx].push(chunk.as_slice());
                }
            }
            let score = |x, y| if x == y { 3 } else { -4 };
            chunks
                .iter()
                .map(|seqs| {
                    let ws = vec![1.; seqs.len()];
                    poa_hmm::POA::generate_w_param(seqs, &ws, -6, -6, &score)
                })
                .collect()
        })
        .collect();
    let config = &poa_hmm::DEFAULT_CONFIG;
    for pos in 0..chain_len {
        for (idx, read) in data.iter().enumerate() {
            let &(_, ref unit) = read.iter().filter(|&&(p, _)| p == pos).nth(0).unwrap();
            let lks: Vec<_> = models
                .iter()
                .map(|ms| ms[pos].forward(unit, config))
                .map(|lk| format!("{}", lk))
                .collect();
            debug!("DUMP\t{}\t{}\t{}", idx, pos, lks.join("\t"));
        }
    }
}
