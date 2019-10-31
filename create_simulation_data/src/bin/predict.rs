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
use bio_utils::fasta::Record;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
const K: usize = 6;
const ORIGINAL: &'static str = "NC_037304.1";
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = bio_utils::fasta::parse_into_vec(&args[1])?;
    debug!("{} reads in total.", reads.len());
    let mut rng: StdRng = SeedableRng::seed_from_u64(23124324);
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[3])?;
    let cores: usize = args[4].parse().unwrap();
    {
        debug!(
            "Could be encoded:{}",
            last_tiling::encoding(&reads, &reference, &alignments).len()
        );
    }
    let reads: Vec<_> = reads
        .into_iter()
        .filter(|e| {
            let desc = e.desc().unwrap();
            !desc.starts_with("junk") && !desc.starts_with("random")
        })
        .collect();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cores)
        .build_global()
        .unwrap();
    debug!("{} alignments in total", alignments.len());
    let len = reference.get_by_id(0).unwrap().len();
    let (training, testset): (Vec<_>, Vec<_>) = reads.into_iter().partition(|r| {
        let desc: Vec<_> = r.desc().unwrap().split(',').collect();
        let forward: bool = desc[1].starts_with('+');
        let region = desc[2].split_whitespace().nth(0).unwrap();
        let start: usize = if forward {
            region.split('-').nth(0).unwrap().parse().unwrap()
        } else {
            len - region.split('-').nth(1).unwrap().parse::<usize>().unwrap()
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
            let is_original = e.desc()?.split(',').nth(0)? == ORIGINAL;
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
            .map(|e| e.unit)
            .min()
            .unwrap_or(0)
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
    Some(predicts)
}

type Seq = Vec<u8>;
// UnitID -> Sequences -> (color, seq)
fn construct_predictor(
    training: Vec<Record>,
    alignments: &[LastTAB],
    contig: &Contigs,
) -> Option<Vec<(Vec<Seq>, Vec<Seq>)>> {
    let (original, mutant): (Vec<_>, Vec<_>) = training
        .into_iter()
        .partition(|read| read.desc().unwrap().split(',').nth(0).unwrap() == ORIGINAL);
    let original = last_tiling::encoding(&original, contig, &alignments);
    let mutant = last_tiling::encoding(&mutant, contig, &alignments);
    let mut chunks = vec![(vec![], vec![]); contig.get_last_unit(0)? as usize + 1];
    use last_tiling::revcmp;
    for read in original {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                revcmp(unit.bases.as_bytes())
            };
            chunks[u].0.push(seq);
        }
    }
    for read in mutant {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                revcmp(unit.bases.as_bytes())
            };
            chunks[u].1.push(seq);
        }
    }
    Some(chunks)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}

fn make_prediction(chunks: &[(Vec<Seq>, Vec<Seq>)], read: &last_tiling::EncodedRead) -> u8 {
    let predicts: Vec<(f64, f64)> = read
        .seq()
        .iter()
        .filter_map(|e| e.encode())
        .map(|e| {
            let u = e.unit as usize;
            let (ref original, ref mutant) = &chunks[u];
            let query = if e.is_forward {
                e.bases.as_bytes().to_vec()
            } else {
                last_tiling::revcmp(e.bases.as_bytes())
            };
            let o = unit_prediction(original, &query);
            let m = unit_prediction(mutant, &query);
            as_weight(o, m)
        })
        .collect();
    // Sum ln2 - H(p)
    let entropy:Vec<_> = predicts
        .iter()
        .map(|(o, m)| (2f64).ln() + o * o.ln() + m * m.ln())
        .collect();
    let sum = entropy.iter().sum::<f64>();
    let (p1, p2) = predicts
        .iter()
        .zip(entropy.iter())
        .map(|((o, m), w)| (o * w, m * w))
        .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
    let (p1, p2) = (p1 / sum, p2 / sum);
    assert!((p1 + p2 - 1.).abs() < 0.0001);
    if p1 < p2 {
        1
    } else {
        0
    }
}

fn unit_prediction(units: &[Seq], query: &[u8]) -> f64 {
    let hmm = DBGHMM::new(&units, K);
    hmm.forward(&query, &DEFAULT_CONFIG)
}
