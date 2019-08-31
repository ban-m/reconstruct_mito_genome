extern crate handmade_bloom_filter;
extern crate kmer_tiling;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate rand;
extern crate serde;
extern crate serde_json;
use bio::io;
use handmade_bloom_filter::{UpperBoundBFFactory, HUGE_MODULO};
use rand::{seq::SliceRandom, thread_rng, Rng};
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Write};

#[derive(Deserialize, Serialize, Debug)]
struct RoundSearch {
    kmer: String,
    k: usize,
    dists: Vec<DistCount>,
}

#[derive(Deserialize, Serialize, Debug)]
struct DistCount {
    dist: u16,
    count: Vec<u16>,
}

fn main() -> std::io::Result<()> {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();
    let args: Vec<_> = std::env::args().collect();
    let t: usize = args[2].parse().unwrap();
    let k: usize = args[3].parse().unwrap();
    let max_dist: u16 = args[4].parse().unwrap();
    let input: Vec<_> = io::fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| e.seq().to_vec())
        .collect();
    info!("Collect reads:{}", input.len());
    let bf = UpperBoundBFFactory::default()
        .k(k)
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    info!("{}", bf);
    let mut rng = thread_rng();
    let result: Vec<_> = (0..10)
        .map(|_| {
            let kmer = loop {
                let line = input.choose(&mut rng).unwrap();
                if line.len() >= k {
                    let start = rng.gen_range(0, line.len() - k);
                    break line[start..start + k].to_vec();
                }
            };
            let dists: Vec<_> = (0..(max_dist + 1))
                .map(|dist| {
                    let p = (0.025f64).powi(dist as i32);
                    let count: Vec<_> =
                        kmer_tiling::enumerate_dist_k_squished_withprob(&kmer, dist, p)
                            .into_iter()
                            .map(|kmer| bf.upper_bound(&kmer).unwrap())
                            .collect();
                    DistCount { dist, count }
                })
                .collect();
            let kmer = String::from_utf8_lossy(&kmer).to_string();
            RoundSearch { dists, kmer, k }
        })
        .collect();
    let mut wtr = BufWriter::new(std::io::stdout());
    let res = serde_json::ser::to_string(&result).unwrap();
    writeln!(&mut wtr, "{}", res)
}
