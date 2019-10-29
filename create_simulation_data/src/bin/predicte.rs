extern crate bio_utils;
extern crate last_tiling;
extern crate rand;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = bio_utils::fasta::parse_into_vec(&args[1])?;
    debug!("{} reads in total.", reads.len());
    let mut rng: StdRng = SeedableRng::seed_from_u64(23124324);
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[3])?;
    let cores: usize = args[4].parse().unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cores)
        .build_global()
        .unwrap();
    debug!("{} alignments in total", alignments.len());
    let len = reference.get_by_id(0).unwrap().len();
    let (training, testset): (Vec<_>, Vec<_>) = reads.into_iter().partition(|r| {
        let desc: Vec<_> = r.desc().unwrap().split(',').collect();
        let forward: bool = desc[1].starts_with('-');
        let start: usize = if forward {
            desc[2].split('-').nth(0).unwrap().parse().unwrap()
        } else {
            len - desc[2].split('-').nth(1).unwrap().parse::<usize>().unwrap()
        };
        if start < len / 2 {
            rng.gen_bool(0.5)
        } else {
            false
        }
    });
    debug!(
        "{} training reads and {} test reads dataset.",
        training.len(),
        testset.len()
    );
    let result = predict(training, testset, alignments, &reference).unwrap();
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    writeln!(&mut wtr, "ReadID\tAnswer\tPredict\tDistance\tLength")?;
    for (readid, answer, predict, distance, length) in result {
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}",
            readid, answer, predict, distance, length
        )?;
    }
    Ok(())
}
use bio_utils::fasta::Record;
use last_tiling::Contigs;
use last_tiling::LastTAB;
fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alingmnets: Vec<LastTAB>,
    contig: &Contigs,
) -> Option<Vec<(String, u8, u8, usize, usize)>> {
    None
}
