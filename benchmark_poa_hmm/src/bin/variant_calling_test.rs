extern crate edlib_sys;
extern crate env_logger;
extern crate last_decompose;
extern crate log;
extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let seed: u64 = std::env::args()
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap_or(132980);
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
    let chain_len = 40;
    let template: Vec<Vec<u8>> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let num_cluster = 2;
    let total_coverage = 300;
    let ws = vec![0.1, 0.9];
    let p = gen_sample::Profile {
        sub: 0.001 / 6.,
        ins: 0.001 / 6.,
        del: 0.001 / 6.,
    };
    let templates: Vec<Vec<_>> = (0..num_cluster)
        .map(|_| {
            template
                .iter()
                .map(|t| gen_sample::introduce_randomness(t, &mut rng, &p))
                .collect()
        })
        .collect();
    let dists: Vec<_> = (0..chain_len)
        .map(|pos| {
            let mut res = 0;
            for i in 0..num_cluster {
                for j in (i + 1)..num_cluster {
                    let t1 = &templates[i][pos];
                    let t2 = &templates[j][pos];
                    res += edlib_sys::global_dist(t1, t2);
                }
            }
            res
        })
        .collect();
    let coverages: Vec<_> = ws
        .iter()
        .map(|r| (r * total_coverage as f64).floor() as usize)
        .collect();
    let mut gen = |ts: &Vec<Vec<u8>>| {
        ts.iter()
            .map(|t| gen_sample::introduce_randomness(&t, &mut rng, &gen_sample::PROFILE))
            .collect::<Vec<_>>()
    };
    let dataset: Vec<_> = coverages
        .iter()
        .enumerate()
        .flat_map(|(idx, &cov)| (0..cov).map(|_| gen(&templates[idx])).collect::<Vec<_>>())
        .collect();
    let score = |x, y| if x == y { 3 } else { -4 };
    let models: Vec<Vec<_>> = {
        let mut chunks = vec![vec![]; chain_len];
        for seq in dataset.iter() {
            for (idx, chunk) in seq.iter().enumerate() {
                chunks[idx].push(chunk.as_slice());
            }
        }
        let weights_of_reads = {
            let cluster: Vec<usize> = (0..num_cluster).collect();
            let mut weights = vec![];
            for _ in 0..total_coverage {
                let mut ws = vec![0.; num_cluster];
                use rand::seq::SliceRandom;
                let ans = cluster.choose(&mut rng).unwrap();
                ws[*ans] = 1.;
                weights.push(ws);
            }
            weights
        };
        let mut weights = vec![vec![]; num_cluster];
        for weights_each_cluster in weights_of_reads {
            for (cluster, &w) in weights_each_cluster.iter().enumerate() {
                weights[cluster].push(w);
            }
        }
        weights
            .iter()
            .map(|ws| {
                chunks
                    .iter()
                    .map(|cs| POA::generate(cs, ws, (-6, -6, &score)))
                    .collect()
            })
            .collect()
    };
    use last_decompose::variant_calling;
    {
        let dataset: Vec<Vec<_>> = dataset
            .iter()
            .map(|read| read.iter().cloned().enumerate().collect())
            .collect();
        let (variants, _) =
            variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, &ws, true);
        for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
            println!("BEFORE\t{}\t{}\t{:.4}\tCentrize", idx, d, w.abs());
        }
        let (variants, _) =
            variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, &ws, false);
        for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
            println!("BEFORE\t{}\t{}\t{:.4}\tNormal", idx, d, w.abs());
        }
    }
    let models: Vec<Vec<_>> = {
        let mut chunks = vec![vec![]; chain_len];
        for seq in dataset.iter() {
            for (idx, chunk) in seq.iter().enumerate() {
                chunks[idx].push(chunk.as_slice());
            }
        }
        let weights_of_reads: Vec<_> = coverages
            .iter()
            .enumerate()
            .flat_map(|(idx, &cov)| {
                let mut weight = vec![0.; num_cluster];
                weight[idx] = 1.;
                vec![weight.clone(); cov]
            })
            .collect();
        let mut weights = vec![vec![]; num_cluster];
        for weights_each_cluster in weights_of_reads {
            for (cluster, &w) in weights_each_cluster.iter().enumerate() {
                weights[cluster].push(w);
            }
        }
        weights
            .iter()
            .map(|ws| {
                chunks
                    .iter()
                    .map(|cs| POA::generate(cs, ws, (-6, -6, &score)))
                    .collect()
            })
            .collect()
    };
    let dataset: Vec<Vec<_>> = dataset
        .into_iter()
        .map(|read| read.into_iter().enumerate().collect())
        .collect();
    let (variants, _) =
        variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, &ws, true);
    for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
        println!("AFTER\t{}\t{}\t{:.4}\tCentrize", idx, d, w.abs());
    }
    let (variants, _) =
        variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, &ws, false);
    for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
        println!("AFTER\t{}\t{}\t{:.4}\tNormal", idx, d, w.abs());
    }
}
