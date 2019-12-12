extern crate bio_utils;
extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate last_tiling;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta::Record;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
const K: usize = 6;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    eprintln!("{:?}", args);
    let reads: Vec<_> = bio_utils::fasta::parse_into_vec(&args[1])?
        .into_iter()
        .filter(|r| r.desc().unwrap().contains("sample1"))
        .collect();
    debug!("{} reads in total.", reads.len());
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[3])?;
    debug!("{} alignments in total", alignments.len());
    let result = coverage_and_likelihood(&reads, alignments, &reference);
    debug!("Dump results");
    for (coverage, likelihood) in result {
        println!("{}\t{}", coverage, likelihood);
    }
    Ok(())
}

fn coverage_and_likelihood(
    reads: &[Record],
    alignments: Vec<LastTAB>,
    contig: &Contigs,
) -> Vec<(f64, f64)> {
    let data = setup(reads, alignments, contig);
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(reads.len() as u64);
    (0..)
        .map(|e| e as f64 * 0.05)
        .take_while(|&p| p < 1.0)
        .flat_map(|pick_prob| {
            let rep_num = pick_prob.recip() as usize;
            (0..rep_num)
                .flat_map(|_| {
                    let (train, test): (Vec<_>, Vec<_>) =
                        data.iter().cloned().partition(|_| rng.gen_bool(pick_prob));
                    cov_and_lk(&train, &test, contig)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn setup(reads: &[Record], alns: Vec<LastTAB>, contig: &Contigs) -> Vec<last_decompose::ERead> {
    let data: Vec<_> = last_tiling::encoding(reads, contig, &alns)
        .into_iter()
        .map(last_decompose::ERead::new)
        .collect();
    data
}

fn cov_and_lk(
    train: &[last_decompose::ERead],
    test: &[last_decompose::ERead],
    contig: &Contigs,
) -> Vec<(f64, f64)> {
    let weights: Vec<_> = train.iter().map(|_| vec![1.]).collect();
    let len: Vec<_> = (0..contig.get_num_of_contigs())
        .map(|e| contig.get_last_unit(e as u16).unwrap() as usize + 1)
        .collect();
    let models = last_decompose::construct_with_weights(train, &weights, K, &len, 0);
    test.par_iter()
        .flat_map(|read| {
            read.seq().par_iter().map(|chunk| {
                let m = &models[chunk.contig as usize][chunk.unit as usize];
                let lk = m.forward(chunk.bases(), &dbg_hmm::PACBIO_CONFIG);
                let w = m.weight();
                (w, lk)
            })
        })
        .collect()
}
