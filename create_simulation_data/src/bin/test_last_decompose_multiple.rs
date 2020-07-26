#[macro_use]
extern crate log;
const LIMIT: u64 = 3600;
use last_decompose::poa_clustering::gibbs_sampling;
use poa_hmm::gen_sample;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    env_logger::init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let (test_num, coverage, probs, clusters, seed, errors) = if args.len() > 4 {
        let tn = args[1].parse::<usize>().unwrap();
        let cov = args[2].parse::<usize>().unwrap();
        let seed = args[3].parse::<u64>().unwrap();
        let errors = args[4].parse::<f64>().unwrap();
        let prob: Vec<_> = args[5..]
            .iter()
            .filter_map(|e| e.parse::<f64>().ok())
            .collect();
        let clusters = prob.len();
        (tn, cov, prob, clusters, seed, errors)
    } else {
        (200, 0, vec![2f64.recip(); 2], 2, 11920981, 0.2)
    };
    let len = 100;
    let chain_len = 90;
    let p = &gen_sample::Profile {
        sub: errors / 6.,
        ins: errors / 6.,
        del: errors / 6.,
    };
    use std::time::Instant;
    let s = Instant::now();
    let (hmm, dists) = benchmark(
        seed, p, coverage, test_num, chain_len, len, &probs, clusters,
    );
    println!("TestNum:{}\tLabeled:{}", test_num, coverage);
    debug!("Elapsed\t{}\t{}", (Instant::now() - s).as_secs(), test_num);
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
    use create_simulation_data::generate_mul_data;
    let profile = &gen_sample::PROFILE;
    let (dataset, label, answer, _border) =
        generate_mul_data(&templates, coverage, test_num, &mut rng, probs, profile);
    let c = &poa_hmm::DEFAULT_CONFIG;
    let data: Vec<Vec<_>> = dataset
        .clone()
        .into_iter()
        .map(|read| read.into_iter().enumerate().collect())
        .collect();
    let forbidden = vec![vec![]; data.len()];
    use last_decompose::poa_clustering::DEFAULT_ALN;
    let coverage = data.iter().map(|r| r.len()).sum::<usize>() / chain_len;
    let pred = gibbs_sampling(
        &data,
        &label,
        Some(&answer),
        &forbidden,
        chain_len,
        clusters,
        LIMIT,
        c,
        &DEFAULT_ALN,
        coverage,
    );
    assert_eq!(pred.len(), label.len() + answer.len());
    let mut result = vec![vec![0; clusters]; clusters];
    //let mut result = vec![vec![0; clusters + 1]; clusters + 1];
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
                            if d != 0 {
                                debug!("{}\t{}", idx, d);
                            }
                            d
                        })
                        .sum::<u32>();
                    debug!("{}\t{}\t{}", i, j, dist);
                    dist
                })
                .collect::<Vec<_>>()
        })
        .collect();
    debug!("{}", dataset.len());
    {
        let probs: Vec<_> = probs.iter().map(|e| format!("{:3}", e)).collect();
        debug!("Probs:[{}]", probs.join(","));
    };
    for i in 0..clusters {
        let tot = answer.iter().filter(|&&e| e as usize == i).count();
        debug!("Cluster {}:{}", i, tot);
    }
    for (pred, ans) in pred.into_iter().zip(label.into_iter().chain(answer)) {
        // let pred = match pred {
        //     Some(res) => res as usize,
        //     None => 0,
        // };
        result[pred as usize][ans as usize] += 1;
    }
    (result, dists)
}
