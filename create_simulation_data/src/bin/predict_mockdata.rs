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
use dbg_hmm::*;
use last_tiling::Contigs;
use last_tiling::EncodedRead;
use last_tiling::LastTAB;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
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
    let len = reference.get_by_id(0).unwrap().len();
    let answer: HashMap<_, _> = reads
        .iter()
        .filter_map(|e| {
            let is_original = e.desc()?.contains("sample1");
            Some((e.id().to_string(), is_original))
        })
        .collect();
    let dist: HashMap<_, _> = reads
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
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1893749823);
    let (training, testset) = reads.into_iter().partition(|_| rng.gen_bool(0.5));
    // let (training, testset): (Vec<_>, Vec<_>) =
    //     reads.into_iter().partition(|r| dist[r.id()] < len / 2);
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
    for (readid, predict, length) in result {
        let answer = answer[&readid];
        let dist = dist[&readid];
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}",
            readid, answer, predict, dist, length
        )?;
    }
    Ok(())
}

#[allow(dead_code)]
fn before_half(r: &Record, len: usize) -> bool {
    let desc: &str = r.desc().unwrap();
    let desc: Vec<_> = desc.split_whitespace().nth(0).unwrap().split(',').collect();
    if desc.len() < 2 {
        debug!("{:?}", r.desc());
    }
    let forward: bool = desc[1].starts_with('+');
    let start: usize = if forward {
        desc[2].split('-').nth(0).unwrap().parse().unwrap()
    } else {
        let l = desc[2].split('-').nth(1).unwrap().parse::<usize>().unwrap();
        len - len.min(l)
    };
    start < len / 2
}

fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
) -> Option<Vec<(String, u8, usize)>> {
    let mut is_original: Vec<_> = training
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
    let mut training: Vec<_> = last_tiling::encoding(&training, contig, &alignments);
    let mut tests: Vec<_> = last_tiling::encoding(&tests, contig, &alignments);
    tests.sort_by_key(|read| {
        read.seq()
            .iter()
            .filter_map(|e| e.encode())
            .map(|e| e.unit)
            .min()
            .unwrap()
    });
    let is_test: std::collections::HashSet<_> = tests.iter().map(|e| e.id().to_string()).collect();
    let len = contig.get_last_unit(0)? as usize + 1;
    while tests.is_empty() {
        assert!(is_original.len() == training.len());
        update(&mut is_original, &mut training, &mut tests, len);
    }
    Some(
        training
            .into_iter()
            .zip(is_original.iter())
            .filter(|(read, predict)| is_test.contains(read.id()))
            .map(|(read, predict)| {
                (
                    read.id().to_string(),
                    if *predict { 0u8 } else { 1u8 },
                    read.recover_raw_sequence().len(),
                )
            })
            .collect(),
    )
}
fn update(
    is_original: &mut Vec<bool>,
    training: &mut Vec<EncodedRead>,
    tests: &mut Vec<EncodedRead>,
    len: usize,
) {
    let tot = training.len() + tests.len();
    let mut chunks = construct_predictor(training, is_original, len);
    debug!("Predictors are constructed:{}", chunks.len());
    for (idx, (os, ms)) in chunks.iter().enumerate() {
        debug!("{}\t{}\t{}", idx, os.len(), ms.len());
    }
    let mut res = vec![];
    while let Some(read) = tests.pop() {
        if let Some(predict) = make_prediction(&chunks, &read) {
            training.push(read);
            is_original.push(predict == 0);
        } else {
            res.push(read);
        }
    }
    tests.append(&mut res);
    assert_eq!(tests.len() + training.len(), tot);
}

fn make_prediction(chunks: &[(Vec<Seq>, Vec<Seq>)], read: &last_tiling::EncodedRead) -> Option<u8> {
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
            let o = unit_predict(&query, &original[..min], K);
            let m = unit_predict(&query, &mutant[..min], K);
            Some((o, m))
        })
        .collect();
    if predicts.is_empty() {
        return None;
    }
    Some(predict_by_naive(&predicts))
}
