#[macro_use]
extern crate log;
// use last_decompose::ERead;
use poa_hmm::gen_sample;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const LIMIT: u64 = 3600;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let args: Vec<_> = std::env::args().collect();
    let (probs, clusters, seed, threads) = {
        let seed = args[1].parse::<u64>().unwrap();
        let threads = args[2].parse::<usize>().unwrap();
        let prob: Vec<_> = args[3..]
            .iter()
            .filter_map(|e| e.parse::<f64>().ok())
            .collect();
        let clusters = prob.len();
        (prob, clusters, seed, threads)
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    debug!("Seed:{}", seed);
    let coverage = 0;
    let chain_len = 90;
    let len = 100;
    // divs, errors, test_nums
    let divs = [1, 3];
    let errors = [0.01, 0.10, 0.15];
    let test_nums = [80, 110, 140];
    let mut params = vec![];
    for &div in divs.iter() {
        for &error in errors.iter() {
            for &test_num in test_nums.iter() {
                params.push((div, error, test_num));
            }
        }
    }
    let result: Vec<_> = params
        .into_par_iter()
        .map(|(div, error, test_num)| {
            let seed = seed + (test_num + div) as u64;
            let seed = seed * (error * 1000f64).floor() as u64;
            use std::time::Instant;
            let s = Instant::now();
            let (hmm, dists) = benchmark(
                seed, div, error, coverage, test_num, chain_len, len, &probs, clusters,
            );
            println!("TestNum:{}\tLabeled:{}", test_num, coverage);
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
            line += &format!("\t{:.3}", error);
            line += &format!("\t{}", div);
            line
        })
        .collect();
    for line in result {
        println!("{}", line);
    }
}

fn benchmark(
    seed: u64,
    div: usize,
    error: f64,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    len: usize,
    probs: &[f64],
    clusters: usize,
) -> (Vec<Vec<u32>>, Vec<Vec<u32>>) {
    let seed = 1003437 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let divs = {
        let mut pos = vec![vec![0; chain_len]; clusters];
        for _ in 0..div {
            let c = rng.gen_range(0, clusters);
            let p = rng.gen_range(0, chain_len);
            pos[c][p] += 1;
        }
        pos
    };
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    let templates: Vec<Vec<_>> = divs
        .into_iter()
        .map(|div_num| {
            template
                .iter()
                .zip(div_num)
                .map(|(seq, e)| {
                    let (mut sub, mut ins, mut del) = (0, 0, 0);
                    for _ in 0..e {
                        match rng.gen_range(0, 3) {
                            0 => sub += 1,
                            1 => ins += 1,
                            2 => del += 1,
                            _ => {}
                        }
                    }
                    gen_sample::introduce_errors(seq, &mut rng, sub, del, ins)
                })
                .collect()
        })
        .collect();
    use create_simulation_data::generate_mul_data;
    let profile = gen_sample::PROFILE.norm().mul(error);
    let (dataset, label, answer, _border) =
        generate_mul_data(&templates, coverage, test_num, &mut rng, probs, &profile);
    let c = &poa_hmm::DEFAULT_CONFIG;
    let data: Vec<Vec<_>> = dataset
        .iter()
        .map(|read| read.iter().map(|e| e.as_slice()).enumerate().collect())
        .collect();
    let forbidden = vec![vec![]; data.len()];
    use last_decompose::poa_clustering::{gibbs_sampling, DEFAULT_ALN};

    let id = rng.gen::<u64>() % 100;
    let config = last_decompose::poa_clustering::ClusteringConfig::new(
        chain_len, clusters, LIMIT, coverage, id, true, c,
    );
    let pred = {
        let answer = Some(answer.as_slice());
        gibbs_sampling(&data, &label, answer, &forbidden, &DEFAULT_ALN, config)
    };
    debug!("Index1\tIndex2\tDist");
    let dists: Vec<Vec<_>> = (0..clusters)
        .map(|i| {
            (i + 1..clusters)
                .map(|j| {
                    let dist = templates[i]
                        .iter()
                        .zip(templates[j].iter())
                        .enumerate()
                        .map(|(idx, (t1, t2))| {
                            let d = bio_utils::alignments::edit_dist(t1, t2);
                            debug!("{}:{}", idx, d);
                            d
                        })
                        .sum::<u32>();
                    debug!("{}\t{}\t{}", i, j, dist);
                    dist
                })
                .collect::<Vec<_>>()
        })
        .collect();
    {
        let probs: Vec<_> = probs.iter().map(|e| format!("{:3}", e)).collect();
        debug!("Probs:[{}]", probs.join(","));
    };
    let mut result = vec![vec![0; clusters]; clusters];
    for (pred, ans) in pred.into_iter().zip(answer) {
        // let pred = match pred {
        //     Some(res) => res as usize,
        //     None => 0,
        // };
        result[pred as usize][ans as usize] += 1;
    }
    (result, dists)
}
