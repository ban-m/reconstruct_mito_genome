extern crate handmade_bloom_filter;
extern crate kmer_tiling;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate rand;
use bio::io;
use handmade_bloom_filter::{UpperBoundBFFactory, HUGE_MODULO};
use rand::{seq::SliceRandom, thread_rng, Rng};
use std::io::{BufWriter, Write};
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();
    let args: Vec<_> = std::env::args().collect();
    let t: usize = args[2].parse().unwrap();
    let input: Vec<_> = io::fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| e.seq().to_vec())
        .collect();
    info!("Collect reads:{}", input.len());
    let k = 20;
    let bf = UpperBoundBFFactory::default()
        .k(k)
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    info!("{}", bf);
    let mut rng = thread_rng();
    let mut wtr = BufWriter::new(std::io::stdout());
    //writeln!(&mut wtr, "ID\tIteration\tCount")?;
    let kmer = loop{
        let line = input.choose(&mut rng).unwrap();
        if line.len() >= k {
            let start = rng.gen_range(0, line.len() - k);
            break line[start..start + k].to_vec();
        }
    };
    let ub = bf.upper_bound(&kmer).unwrap();
    writeln!(&mut wtr, "SeedKmer:{}\tCount:{}", String::from_utf8_lossy(&kmer), ub)?;
    writeln!(&mut wtr, "Kmer\tDist\tCount")?;
    for dist in 1..=2 {
        for kmer in kmer_tiling::enumerate_dist_k(&kmer, dist) {
            let ub = bf.upper_bound(&kmer).unwrap();
            let kmer = String::from_utf8_lossy(&kmer);
            writeln!(&mut wtr, "{}\t{}\t{}", kmer , dist, ub)?;
        }
    }
    Ok(())
}
