extern crate handmade_bloom_filter;
extern crate kmer_tiling;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use bio::io;
use handmade_bloom_filter::{UpperBoundBFFactory, HUGE_MODULO};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

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
    let assignment = open_assignment(&args[5])?;
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
    let input: Vec<_> = io::fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| assignment.contains_key(e.id()))
        .collect();
    for rec in input {
        let line = rec.seq();
        let len = line.len();
        let clip = len / 10;
        let result: Vec<_> = line
            .windows(k)
            .skip(clip)
            .take((clip * 3).min(20))
            .map(|kmer| {
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
        let name = format!("./result/{}_{}.json", k, rec.id().replace("/", ""));
        let mut wtr = BufWriter::new(File::create(Path::new(&name))?);
        let res = serde_json::ser::to_string(&result)?;
        writeln!(&mut wtr, "{}", res)?;
    }
    Ok(())
}

fn open_assignment(file: &str) -> std::io::Result<HashMap<String, String>> {
    use std::io::{BufRead, BufReader};
    let mut res = HashMap::new();
    for line in BufReader::new(File::open(&Path::new(file))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        let mut contents = line.split('\t');
        let id = contents.next().unwrap().to_string();
        let chrtype = contents.next().unwrap();
        let chrtype = match chrtype.parse::<u8>() {
            Ok(_) => "genomic".to_string(),
            Err(_) => chrtype.to_string(),
        };
        res.insert(id, chrtype);
    }
    Ok(res)
}
