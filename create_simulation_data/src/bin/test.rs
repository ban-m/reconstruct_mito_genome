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
// use rayon::prelude::*;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (test_num, coverage, _prob) = if args.len() > 3 {
        let tn = args[1].parse::<usize>().unwrap();
        let cov = args[2].parse::<usize>().unwrap();
        let prob = args[3].parse::<f64>().unwrap();
        (tn, cov, prob)
    } else {
        (200, 0, 0.5)
    };
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let seed = 11920981;
    benchmark(seed, p, coverage, test_num, chain_len, k, len);
}

fn benchmark(
    seed: u64,
    p: &gen_sample::Profile,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    k: usize,
    len: usize,
) {
    let seed = 1003437 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let template1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let template2: Vec<_> = template1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let (dataset, _label, _answer, _border) = generate_mul_data(
        &[template1.clone(), template2.clone()],
        coverage,
        test_num,
        &mut rng,
        &[0.5, 0.5],
    );
    let data: Vec<_> = dataset
        .into_iter()
        .enumerate()
        .map(|(idx, e)| {
            let id = format!("{}", idx);
            last_decompose::ERead::new_with_lowseq(e, &id)
        })
        .collect();
    let contigs = vec![chain_len];
    let mut mf = last_decompose::ModelFactory::new(&contigs, &data, k);
    let weights: Vec<_> = (0..data.len()).map(|_| vec![0.5, 0.5]).collect();
    let models = mf.generate_model(&weights, &data, 0);
    for read in &data {
        for u in read.seq.iter() {
            let model = &models[u.contig()][u.unit()];
            let lk = model.forward(u.bases(), &dbg_hmm::DEFAULT_CONFIG);
            debug!("LK\t{:.3}\t{:.3}", model.weight(), lk);
        }
    }
    let mut f = Factory::new();
    let weights = vec![1.; data.len()];
    let chunks = {
        let mut res = vec![vec![]; chain_len];
        for read in &data {
            for chunk in read.seq.iter() {
                res[chunk.unit()].push(chunk.bases());
            }
        }
        res
    };
    let mut buf = vec![];
    let models: Vec<_> = chunks
        .iter()
        .map(|chunks| f.generate_with_weight_prior(chunks, &weights, k, &mut buf))
        .collect();
    for m in &models {
        debug!("Model\t{}", m);
    }
    for read in &data {
        for u in read.seq.iter() {
            let model = &models[u.unit()];
            let lk = model.forward(u.bases(), &DEFAULT_CONFIG);
            debug!("NN\t{:.3}\t{:.3}", model.weight(), lk);
        }
    }
    let i = 3;
    let chunks: Vec<_> = data.iter().map(|read| read.seq[i].bases()).collect();
    let m = f.generate_with_weight_prior(&chunks, &vec![0.5; chunks.len()], k, &mut buf);
    debug!("TT\t{}", m);
    for e in 0..100 {
        let temp = if e % 2 == 0 {
            &template1[i]
        } else {
            &template2[i]
        };
        let q = gen_sample::introduce_randomness(temp, &mut rng, &gen_sample::PROFILE);
        let lk = m.forward(&q, &DEFAULT_CONFIG);
        debug!("TT\t{:.3}\t{:.3}", m.weight(), lk);
    }
}
