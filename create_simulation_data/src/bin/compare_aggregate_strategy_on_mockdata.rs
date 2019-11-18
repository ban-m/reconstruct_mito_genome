extern crate bio_utils;
extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate last_tiling;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta::Record;
use create_simulation_data::*;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::time;
const K: usize = 6;
// const MAX_COVERAGE: usize = 30;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("predict_mockdata=debug"))
        .init();
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = bio_utils::fasta::parse_into_vec(&args[1])?
        .into_iter()
        .filter(|r| r.desc().unwrap().contains("sample"))
        .collect();
    debug!("{} reads in total.", reads.len());
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[3])?;
    debug!("{} alignments in total", alignments.len());
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1893749823);
    let (training, testset): (Vec<_>, Vec<_>) = reads.into_iter().partition(|_| rng.gen_bool(0.5));
    let training: Vec<_> = training.into_iter().filter(|_| rng.gen_bool(0.4)).collect();
    debug!(
        "{} training reads and {} test reads dataset.",
        training.len(),
        testset.len()
    );
    let result = predict(training, testset, alignments, &reference).unwrap();
    debug!("Dump results");
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    writeln!(&mut wtr, "ReadID\tAnswer\tPredict\tLength\tDistance")?;
    for (readid, answer, predict, length, dist) in result {
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}",
            readid, answer, predict, dist, length
        )?;
    }
    Ok(())
}

fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
) -> Option<Vec<(String, u8, u8, usize, usize)>> {
    let answer: HashMap<_, _> = tests
        .iter()
        .filter_map(|e| {
            let is_original = e.desc()?.contains("sample1");
            Some((e.id().to_string(), is_original))
        })
        .collect();
    let len = contig.get_by_id(0)?.len();
    let dist: HashMap<_, _> = tests
        .iter()
        .filter_map(|e| {
            let desc: Vec<_> = e
                .desc()?
                .split_whitespace()
                .nth(0)
                .unwrap()
                .split(',')
                .collect();
            let forward: bool = desc[1].starts_with('+');
            let dist = if forward {
                desc[2].split('-').nth(0)?.parse().ok()?
            } else {
                let l = desc[2].split('-').nth(1)?.parse::<usize>().ok()?;
                len - l.min(len)
            };
            Some((e.id().to_string(), dist))
        })
        .collect();
    debug!("Answer and Distance hashmaps are built");
    let mut tests = last_tiling::encoding(&tests, contig, &alignments);
    tests.sort_by_key(|read| {
        let (min, _max) = read
            .seq()
            .iter()
            .filter_map(|e| e.encode())
            .map(|e| e.unit)
            .fold((std::u16::MAX, std::u16::MIN), |(x, y), u| {
                (x.min(u), y.max(u))
            });
        (min, read.seq().len())
    });
    debug!("Tetst cases are sorted.");
    let mut predicts = vec![];
    let is_original: Vec<_> = training
        .iter()
        .map(|read| {
            read.desc()
                .unwrap()
                .split(',')
                .nth(0)
                .unwrap()
                .contains("sample1")
        })
        .collect();
    let training = last_tiling::encoding(&training, contig, &alignments);
    let len = contig.get_last_unit(0)? as usize + 1;
    let mut chunks = construct_predictor(&training, &is_original, len);
    debug!("Predictors are constructed:{}", chunks.len());
    for (idx, (os, ms)) in chunks.iter().enumerate() {
        debug!("{}\t{}\t{}", idx, os.len(), ms.len());
    }
    for read in tests {
        let id = read.id().to_string();
        let answer = if answer[&id] { 0 } else { 1 };
        let length = read.recover_raw_sequence().len();
        let dist = dist[&id];
        let predict = make_prediction(&chunks, &read);
        debug!("Pred\tAns={}\t{}", predict, answer);
        predicts.push((id, answer, predict, length, dist));
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                last_tiling::revcmp(unit.bases.as_bytes())
            };
            if predict == 0 {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    debug!("Finish predicting");
    Some(predicts)
}

fn make_prediction(chunks: &[(Vec<Seq>, Vec<Seq>)], read: &last_tiling::EncodedRead) -> u8 {
    let start = time::Instant::now();
    let predicts: Vec<(f64, f64)> = read
        .seq()
        .par_iter()
        .filter_map(|e| e.encode())
        .filter_map(|e| {
            let u = e.unit as usize;
            let (ref original, ref mutant) = &chunks[u];
            if original.len() < 10 || mutant.len() < 10 {
                return None;
            }
            let query = if e.is_forward {
                e.bases.as_bytes().to_vec()
            } else {
                last_tiling::revcmp(e.bases.as_bytes())
            };
            let min = original.len().min(mutant.len());
            // let s = time::Instant::now();
            let o = unit_predict(&query, &original[..min], K);
            let m = unit_predict(&query, &mutant[..min], K);
            Some((o, m))
        })
        .collect();
    if predicts.is_empty() {
        debug!("Error! Overall:{:.4}\t{:.4}", 0.5, 0.5);
        return 0;
    }
    let res = predict_by_naive(&predicts);
    // let res = predict_by_sow(&predicts);
    // let res = predict_by_sol(&predicts);
    debug!("Time {}", (time::Instant::now() - start).as_millis());
    res
}
