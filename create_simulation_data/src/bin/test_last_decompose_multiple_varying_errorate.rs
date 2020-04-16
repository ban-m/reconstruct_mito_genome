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
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
const LIMIT: u64 = 3600;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(8)
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
    let chain_len = 90;
    let len = 100;
    let errors = [1, 1, 2, 2, 3, 3, 4, 4];
    let test_nums: Vec<_> = (80..=150).step_by(10).collect();
    for &error in errors.iter() {
        for &test_num in &test_nums {
            let seed = seed + (test_num + error) as u64;
            println!("TestNum:{}\tLabeled:{}", test_num, coverage);
            use std::time::Instant;
            let s = Instant::now();
            let (hmm, dists) = benchmark(
                seed, error, coverage, test_num, chain_len, len, &probs, clusters,
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
            line += &format!("\t{}", error);
            println!("{}", line);
        }
    }
}

fn benchmark(
    seed: u64,
    error: usize,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    len: usize,
    probs: &[f64],
    clusters: usize,
) -> (Vec<Vec<u32>>, Vec<Vec<u32>>) {
    let seed = 1003437 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let errors = {
        let mut pos = vec![vec![0; chain_len]; clusters];
        for _ in 0..error {
            let c = rng.gen_range(0, clusters);
            let p = rng.gen_range(0, chain_len);
            pos[c][p] += 1;
        }
        pos
    };
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    let templates: Vec<Vec<_>> = errors
        .into_iter()
        .map(|error_num| {
            template
                .iter()
                .zip(error_num)
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
                            let d = edlib_sys::global_dist(t1, t2);
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
    let _ok_chunk: Vec<bool> = (0..chain_len)
        .map(|idx| {
            (0..clusters).any(|i| {
                (i + 1..clusters).any(|j| {
                    let d = edlib_sys::global_dist(&templates[i][idx], &templates[j][idx]);
                    d != 0
                })
            })
        })
        .collect();
    let (dataset, label, answer, _border) =
        create_simulation_data::generate_mul_data(&templates, coverage, test_num, &mut rng, probs);
    let c = &poa_hmm::DEFAULT_CONFIG;
    {
        let probs: Vec<_> = probs.iter().map(|e| format!("{:3}", e)).collect();
        debug!("Probs:[{}]", probs.join(","));
    };
    let data: Vec<_> = dataset
        .into_iter()
        .enumerate()
        .map(|(idx, e)| {
            let id = format!("{}", idx);
            let read = ERead::new_with_lowseq(e, &id);
            read
        })
        .collect();
    let forbidden = vec![vec![]; data.len()];
    let em_pred = clustering(&data, (&label, &answer), &forbidden, clusters, LIMIT, c);
    let mut result = vec![vec![0; clusters]; clusters];
    for (pred, ans) in em_pred.into_iter().zip(answer) {
        let pred = match pred {
            Some(res) => res as usize,
            None => 0,
        };
        result[pred][ans as usize] += 1;
    }
    (result, dists)
}
