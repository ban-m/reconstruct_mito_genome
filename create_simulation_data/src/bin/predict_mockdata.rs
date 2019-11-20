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
use last_tiling::EncodedRead;
use last_tiling::LastTAB;
use rand::{Rng,SeedableRng};
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
    let (training, testset): (Vec<_>, Vec<_>) = reads.into_iter().partition(|_| rng.gen_bool(0.2));
    // let (training, testset): (Vec<_>, Vec<_>) = reads
    //     .into_iter()
    //     .partition(|r| dist[r.id()] < len / 2);
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
    let (mut sum, mut ok) = (0, 0);
    for (readid, predict, length) in result {
        let answer = if answer[&readid] { 0 } else { 1 };
        let dist = dist[&readid];
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}",
            readid, answer, predict, dist, length
        )?;
        sum += 1;
        ok += if answer == predict { 1 } else { 0 };
    }
    debug!("Result:{}\t{}", sum, ok);
    Ok(())
}

fn setup(
    training: Vec<Record>,
    tests: Vec<Record>,
    alns: Vec<LastTAB>,
    contig: &Contigs,
) -> (Vec<bool>, Vec<bool>, Vec<EncodedRead>, usize) {
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
        .chain(vec![false; tests.len()])
        .collect();
    let mut training: Vec<_> = last_tiling::encoding(&training, &contig, &alns);
    let mut tests: Vec<_> = last_tiling::encoding(&tests, contig, &alns);
    tests.sort_by_key(|read| {
        read.seq()
            .iter()
            .filter_map(|e| e.encode())
            .map(|e| e.unit)
            .min()
            .unwrap()
    });
    let mut background = vec![false; training.len()];
    for _ in 0..tests.len() / 2 {
        background.push(false);
        background.push(true);
    }
    if background.len() < is_original.len() {
        background.push(false);
    }
    let training_size = training.len();
    training.append(&mut tests);
    (is_original, background, training, training_size)
}

fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
) -> Option<Vec<(String, u8, usize)>> {
    // let is_test: std::collections::HashSet<_> = tests.iter().map(|e| e.id().to_string()).collect();
    let (mut is_original, mut bg_assign, data, border) = setup(training, tests, alignments, contig);
    let len = contig.get_last_unit(0)? as usize + 1;
    let tot = data.len();
    let ori = is_original.iter().filter(|&&e| e).count();
    let muta = tot - ori;
    debug!("Begin!");
    debug!("{}\t{}\t{}", tot, ori, muta);
    while update(&mut is_original, &mut bg_assign, &data, border, len) {
        let tot = data.len();
        let ori = is_original.iter().filter(|&&e| e).count();
        let muta = tot - ori;
        debug!("Updated");
        debug!("{}\t{}\t{}", tot, ori, muta);
    }
    let result: Vec<_> = data[border..]
        .iter()
        .zip(is_original[border..].iter())
        .map(|(read, predict)| {
            let id = read.id().to_string();
            let pred = if *predict { 0u8 } else { 1u8 };
            let len = read.recover_raw_sequence().len();
            (id, pred, len)
        })
        .collect();
    Some(result)
}
fn update(
    is_original: &mut Vec<bool>,
    bg_assign: &mut Vec<bool>,
    data: &[EncodedRead],
    border: usize,
    len: usize,
) -> bool {
    use std::time::Instant;
    let s = Instant::now();
    // let mut chunks = construct_predictors(&data, is_original, bg_assign, len);
    let mut chunks = construct_predictors_simple(&data, is_original, len);
    debug!("Predictors are constructed:{}", chunks.len());
    // for (idx, (os, ms)) in chunks.iter().enumerate() {
    //     if idx % 4 == 0 {
    //         debug!("{}\t{}\t{}", idx, os.len(), ms.len());
    //     }
    // }
    let mut state = false;
    let mut is_updated = false;
    is_original[border..]
        .iter_mut()
        .zip(data[border..].iter())
        .zip(bg_assign[border..].iter_mut())
        .filter(|&((&mut is_ori, _), _)| !is_ori)
        .for_each(|((is_ori, read), bg)| {
            // let predict = make_prediction(&chunks, &read, *bg);
            let predict = make_prediction_simple(&chunks, &read);
            is_updated |= predict != if *is_ori { 0 } else { 1 };
            if predict == 0 {
                //merge(&mut chunks, read);
                // merge_simple(&mut chunks, read);
                *is_ori = true;
            } else {
                *is_ori = false;
                *bg = state;
                state = !state;
            }
        });
    debug!("{:?}", Instant::now() - s);
    is_updated
}

fn make_prediction_simple(chunks: &[(Vec<Seq>, Vec<Seq>)], read: &last_tiling::EncodedRead) -> u8 {
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
    //debug!("{}\t{}", read.id(), predicts.len());
    if predicts.len() < 10 {
        return 1;
    }
    predict_by_naive(&predicts)
}

fn construct_predictors_simple(
    data: &[EncodedRead],
    is_original: &[bool],
    len: usize,
) -> Vec<(Vec<Seq>, Vec<Seq>)> {
    let mut chunks = vec![(vec![], vec![]); len];
    use last_tiling::revcmp;
    for (&is_original, read) in is_original.iter().zip(data.iter()) {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                revcmp(unit.bases.as_bytes())
            };
            if is_original {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    chunks
}

