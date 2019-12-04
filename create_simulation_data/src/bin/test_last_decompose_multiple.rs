extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::gen_sample;
use last_decompose::{clustering, likelihood_of_assignments, ERead};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (test_num, coverage, probs, clusters) = if args.len() > 3 {
        let tn = args[1].parse::<usize>().unwrap();
        let cov = args[2].parse::<usize>().unwrap();
        let prob: Vec<_> = args[3..]
            .iter()
            .filter_map(|e| e.parse::<f64>().ok())
            .collect();
        let clusters = prob.len();
        (tn, cov, prob, clusters)
    } else {
        (200, 0, vec![2f64.recip(); 2], 2)
    };
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    use std::time::Instant;
    let seed = 11920981;
    let s = Instant::now();
    let (hmm, dists) = benchmark(
        seed, p, coverage, test_num, chain_len, k, len, &probs, clusters,
    );
    eprintln!("{:?}", Instant::now() - s);
    info!("Cov:{}", coverage);
    for (idx, preds) in hmm.into_iter().enumerate() {
        let tp = preds[idx];
        let tot = preds.iter().sum::<u32>();
        eprint!("Predicted as {}:", idx);
        for ans in preds {
            eprint!("{}\t", ans);
        }
        eprintln!("Total:{:.4}", tp as f64 / tot as f64);
    }
    for (idx, ds) in dists.into_iter().enumerate() {
        eprint!("Distance from {}:", idx);
        for d in ds {
            eprint!("{}\t", d);
        }
        eprintln!();
    }
}

fn benchmark(
    seed: u64,
    p: &gen_sample::Profile,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    k: usize,
    len: usize,
    probs: &[f64],
    clusters: usize,
) -> (Vec<Vec<u32>>, Vec<Vec<u32>>) {
    let seed = 1003437 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut templates = vec![];
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    for _ in 0..clusters {
        let seq: Vec<_> = template
            .iter()
            .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
            .collect();
        templates.push(seq);
    }
    debug!("Index1\tIndex2\tDist");
    let dists: Vec<Vec<_>> = (0..clusters)
        .map(|i| {
            (i + 1..clusters)
                .map(|j| {
                    let dist = templates[i]
                        .iter()
                        .zip(templates[j].iter())
                        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
                        .sum::<u32>();
                    debug!("{}\t{}\t{}", i, j, dist);
                    dist
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let (dataset, label, answer, _border) =
        create_simulation_data::generate_mul_data(&templates, coverage, test_num, &mut rng, probs);
    let contigs = vec![chain_len];
    let data: Vec<_> = dataset
        .into_iter()
        .enumerate()
        .map(|(idx, e)| {
            let id = format!("{}", idx);
            ERead::new_with_lowseq(e, &id)
        })
        .collect();
    {
        let answer: Vec<_> = label.iter().chain(answer.iter()).copied().collect();
        let objlk = likelihood_of_assignments(&data, &answer, k, clusters, &contigs);
        debug!("ObjLK:{}", objlk);
    }
    let forbidden = vec![vec![]; data.len()];
    let em_pred = clustering(&data, &label, &forbidden, k, clusters, &contigs);
    let mut result = vec![vec![0; clusters]; clusters];
    for (ans, pred) in em_pred.into_iter().zip(answer) {
        result[pred as usize][ans as usize] += 1;
    }
    (result, dists)
}
