extern crate bio_utils;
extern crate dbg_hmm;
extern crate last_tiling;
extern crate rand;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use dbg_hmm::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
// use rayon::prelude::*;
use std::io::{BufWriter, Write};
use bio_utils::fasta::Record;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use std::collections::HashMap;
const K: usize = 6;
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
fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
) -> Option<Vec<(String, u8, u8, usize, usize)>> {
    let answer: HashMap<_, _> = tests
        .iter()
        .filter_map(|e| {
            let is_original = e.desc()?.split(',').nth(0)? == "original";
            Some((e.id().to_string(), is_original))
        })
        .collect();
    let len = contig.get_by_id(0)?.len();
    let dist: HashMap<_, _> = tests
        .iter()
        .filter_map(|e| {
            let desc: Vec<_> = e.desc()?.split(',').collect();
            let forward: bool = desc[1].starts_with('-');
            let dist = if forward {
                desc[2].split('-').nth(0)?.parse().ok()?
            } else {
                len - desc[2].split('-').nth(1)?.parse::<usize>().ok()?
            };
            Some((e.id().to_string(), dist))
        })
        .collect();
    let mut tests = last_tiling::encoding(&tests, contig, &alignments);
    tests.sort_by_key(|read| {
        read.seq()
            .iter()
            .filter_map(|e| e.encode())
            .map(|e| (e.contig, e.unit))
            .fold((100, 2000), |c, t| if c < t { c } else { t })
    });
    let mut predicts = vec![];
    let mut chunks = construct_predictor(training, &alignments, contig)?;
    for read in tests {
        let id = read.id().to_string();
        let answer = if answer[&id] { 0 } else { 1 };
        let length = read.recover_raw_sequence().len();
        let dist = dist[&id];
        let predict = make_prediction(&chunks, &read);
        predicts.push((id, answer, predict, length, dist));
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            if unit.is_forward {
                chunks[u].push((predict, unit.bases.as_bytes().to_vec()));
            } else {
                chunks[u].push((predict, last_tiling::revcmp(unit.bases.as_bytes())));
            }
        }
    }
    Some(predicts)
}

// UnitID -> Sequences -> (color, seq)
fn construct_predictor(
    training: Vec<Record>,
    alignments: &[LastTAB],
    contig: &Contigs,
) -> Option<Vec<Vec<(u8, Vec<u8>)>>> {
    let (original, mutant): (Vec<_>, Vec<_>) = training
        .into_iter()
        .partition(|read| read.desc().unwrap().split(',').nth(0).unwrap() == "original");
    let original = last_tiling::encoding(&original, contig, &alignments);
    let mutant = last_tiling::encoding(&mutant, contig, &alignments);
    let mut chunks = vec![vec![]; contig.get_last_unit(0)? as usize + 1];
    use last_tiling::revcmp;
    for (class, read) in original
        .into_iter()
        .map(|e| (0, e))
        .chain(mutant.into_iter().map(|e| (1, e)))
    {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            if unit.is_forward {
                chunks[u].push((class, unit.bases.as_bytes().to_vec()));
            } else {
                chunks[u].push((class, revcmp(unit.bases.as_bytes())));
            }
        }
    }
    Some(chunks)
}

fn make_prediction(chunks: &[Vec<(u8, Vec<u8>)>], read: &last_tiling::EncodedRead) -> u8 {
    let predicts: Vec<_> = read
        .seq()
        .iter()
        .filter_map(|e| e.encode())
        .map(|e| {
            let u = e.unit as usize;
            let units = &chunks[u];
            if e.is_forward {
                unit_prediction(units, last_tiling::revcmp(e.bases.as_bytes()))
            } else {
                unit_prediction(units, e.bases.as_bytes().to_vec())
            }
        })
        .collect();
    // Sum ln2 - H(p)
    let sum_of_entropy = predicts
        .iter()
        .map(|ps| (2f64).ln() - ps.iter().map(|e| -e * e.ln()).fold(0., |x, y| x + y))
        .fold(0., |x, y| x + y);
    let sum = predicts
        .iter()
        .map(|ps| {
            let weight: f64 = (2f64).ln() - ps.iter().map(|e| -e * e.ln()).fold(0., |x, y| x + y);
            ps.iter().map(|e| e * weight).collect::<Vec<_>>()
        })
        .fold(vec![0., 0.], |acc, ps| {
            acc.into_iter()
                .zip(ps.into_iter())
                .map(|(x, y)| x + y)
                .collect::<Vec<_>>()
        });
    let (arg, _max) = sum
        .into_iter()
        .map(|e| e / sum_of_entropy)
        .enumerate()
        .fold(
            (0, -1.),
            |(arg, max), (idx, prob)| if max < prob { (idx, prob) } else { (arg, max) },
        );
    arg as u8
}

fn unit_prediction(units: &[(u8, Vec<u8>)], query: Vec<u8>) -> Vec<f64> {
    let mut colors: Vec<_> = units.iter().map(|e| e.0).collect();
    colors.sort();
    colors.dedup();
    let probs: Vec<_> = colors
        .iter()
        .map(|color| {
            let train: Vec<&[u8]> = units
                .iter()
                .filter_map(|(i, x)| if i == color { Some(x.as_slice()) } else { None })
                .collect();
            let hmm = DBGHMM::new_from_ref(&train, K);
            hmm.forward(&query, &DEFAULT_CONFIG)
        })
        .collect();
    normalize(probs)
}

fn normalize(probs: Vec<f64>) -> Vec<f64> {
    // Log sum product
    let max = probs
        .iter()
        .fold(std::f64::MIN, |x, &y| if x < y { y } else { x });
    let log_denominator = probs
        .iter()
        .map(|&e| e - max)
        .map(|x| x.exp())
        .fold(0., |x, y| x + y)
        .ln()
        + max;
    probs
        .into_iter()
        .map(|x| x - log_denominator)
        .map(|x| x.exp())
        .collect()
}

