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
use last_decompose::{clustering, ERead};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (probs, clusters, seed) = {
        let seed = args[1].parse::<u64>().unwrap();
        let prob: Vec<_> = args[2..]
            .iter()
            .filter_map(|e| e.parse::<f64>().ok())
            .collect();
        let clusters = prob.len();
        (prob, clusters, seed)
    };
    debug!("Seed:{}", seed);
    let coverage = 0;
    let chain_len = 40;
    let k = 6;
    let len = 150;
    let prob: Vec<_> = (1..=20).step_by(2).collect();
    let test_nums: Vec<_> = (50..201).step_by(30).collect();
    for p in prob {
        for &test_num in &test_nums {
            let p = p as f64 / 6000.;
            println!("ErrorRate:{:3}", p * 6.);
            let p = &gen_sample::Profile {
                sub: p,
                ins: p,
                del: p,
            };
            use std::time::Instant;
            println!("TestNum:{}\tLabeled:{}", test_num, coverage);
            let s = Instant::now();
            let (hmm, dists) = benchmark(
                seed, p, coverage, test_num, chain_len, k, len, &probs, clusters,
            );
            debug!("Elapsed {:?}", Instant::now() - s);
            let mut line = "RESULT".to_string();
            let num_of_reads: Vec<_> = (0..clusters)
                .map(|cl| hmm.iter().map(|pred| pred[cl]).sum::<u32>())
                .collect();
            for nor in &num_of_reads {
                line += &format!("\t{}", nor);
            }
            for (idx, preds) in hmm.into_iter().enumerate() {
                let tp = preds[idx];
                let tot = num_of_reads[idx];
                print!("Predicted as {}:", idx);
                for ans in preds {
                    print!("{}\t", ans);
                    line += &format!("\t{}", ans);
                }
                println!("Total:{:.4}", tp as f64 / tot as f64);
            }
            for (idx, ds) in dists.into_iter().enumerate() {
                print!("Distance from {}:", idx);
                for d in ds {
                    line += &format!("\t{}", d);
                    print!("{}\t", d);
                }
                println!();
            }
            line += &format!("\t{}", test_num);
            println!("{}", line);
        }
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
    let c = &dbg_hmm::DEFAULT_CONFIG;
    {
        let probs: Vec<_> = probs.iter().map(|e| format!("{:3}", e)).collect();
        debug!("Probs:[{}]", probs.join(","));
    };
    let data: Vec<_> = dataset
        .into_iter()
        .enumerate()
        .map(|(idx, e)| {
            let id = format!("{}", idx);
            ERead::new_with_lowseq(e, &id)
        })
        .collect();
    let forbidden = vec![vec![]; data.len()];
    let em_pred = clustering(&data, &label, &forbidden, k, clusters, &contigs, &answer, c);
    let mut result = vec![vec![0; clusters]; clusters];
    for (pred, ans) in em_pred.into_iter().zip(answer) {
        result[pred as usize][ans as usize] += 1;
    }
    (result, dists)
}
