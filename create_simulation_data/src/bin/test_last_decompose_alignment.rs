extern crate bio;
extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rayon;
use dbg_hmm::gen_sample;
use last_decompose::clustering_via_alignment;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (test_num, coverage, probs, clusters, seed) = if args.len() > 3 {
        let tn = args[1].parse::<usize>().unwrap();
        let cov = args[2].parse::<usize>().unwrap();
        let seed = args[3].parse::<u64>().unwrap();
        let prob: Vec<_> = args[4..]
            .iter()
            .filter_map(|e| e.parse::<f64>().ok())
            .collect();
        let clusters = prob.len();
        (tn, cov, prob, clusters, seed)
    } else {
        (200, 0, vec![2f64.recip(); 2], 2, 11920981)
    };
    let chain_len = 20;
    let len = 150;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    use std::time::Instant;
    println!("TestNum:{}\tLabeled:{}", test_num, coverage);
    let s = Instant::now();
    let (hmm, dists) = benchmark(
        seed, p, coverage, test_num, chain_len, len, &probs, clusters,
    );
    debug!("Elapsed {:?}", Instant::now() - s);
    for (idx, preds) in hmm.into_iter().enumerate() {
        let tp = preds[idx];
        let tot = preds.iter().sum::<u32>();
        print!("Predicted as {}:", idx);
        for ans in preds {
            print!("{}\t", ans);
        }
        println!("Total:{:.4}", tp as f64 / tot as f64);
    }
    for (idx, ds) in dists.into_iter().enumerate() {
        print!("Distance from {}:", idx);
        for d in ds {
            print!("{}\t", d);
        }
        println!();
    }
}

fn benchmark(
    seed: u64,
    p: &gen_sample::Profile,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
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
    let reads: Vec<Vec<u8>> = dataset
        .into_iter()
        .map(|e| e.into_iter().flat_map(|e| e).collect::<Vec<u8>>())
        .collect();
    let forbidden = vec![vec![]; reads.len()];
    let read_index: Vec<_> = (0..reads.len()).collect();
    use std::collections::HashMap;
    let sims: Vec<HashMap<usize, i64>> = reads
        .par_iter()
        .map(|r| {
            reads
                .iter()
                .enumerate()
                .map(|(i, q)| (i, alignment(r, q)))
                .collect()
        })
        .collect();
    let aln_pred = clustering_via_alignment(&read_index, &label, &forbidden, &sims, clusters);
    debug!("Finish alignment.");
    let mut result = vec![vec![0; clusters]; clusters];
    for i in 0..clusters {
        let tot = answer.iter().filter(|&&e| e as usize == i).count();
        debug!("Cluster {}:{}", i, tot);
    }
    for (pred, ans) in aln_pred.into_iter().zip(answer) {
        result[pred as usize][ans as usize] += 1;
    }
    (result, dists)
}

fn alignment(r: &[u8], q: &[u8]) -> i64 {
    use bio::alignment::pairwise::Aligner;
    let mut a = Aligner::new(-3, -1, |a, b| if a == b { 1 } else { -1 });
    a.global(r, q).score as i64
}
