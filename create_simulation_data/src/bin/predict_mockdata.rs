extern crate bio_utils;
extern crate dbg_hmm;
extern crate last_tiling;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use bio_utils::fasta::Record;
use dbg_hmm::*;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::time;
const K: usize = 6;
fn main() -> std::io::Result<()>{
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
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
    let (training, testset): (Vec<_>, Vec<_>) =
        reads.into_iter().partition(|r| before_half(r, len));
    debug!(
        "{} training reads and {} test reads dataset.",
        training.len(),
        testset.len()
    );
    let result = predict(training, testset, alignments, &reference).unwrap();
    debug!("Dump results");
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

fn before_half(r: &Record, len: usize) -> bool {
    let desc: &str = r.desc().unwrap();
    let desc: Vec<_> = desc.split_whitespace().nth(0).unwrap().split(',').collect();
    let forward: bool = desc[1].starts_with('+');
    let start: usize = if forward {
        desc[2].split('-').nth(0).unwrap().parse().unwrap()
    } else {
        len - desc[2].split('-').nth(1).unwrap().parse::<usize>().unwrap()
    };
    start < len / 2
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
            // debug!("{:?}",desc[2].split('-').nth(1));
            let dist = if forward {
                desc[2].split('-').nth(0)?.parse().ok()?
            } else {
                len - desc[2].split('-').nth(1)?.parse::<usize>().ok()?
            };
            Some((e.id().to_string(), dist))
        })
        .collect();
    assert!(
        answer.len() == dist.len(),
        "{}\t{}",
        answer.len(),
        dist.len()
    );
    debug!("Answer and Distance hashmaps are built");
    let mut tests = last_tiling::encoding(&tests, contig, &alignments);
    tests.sort_by_key(|read| {
        read.seq()
            .iter()
            .filter_map(|e| e.encode())
            .map(|e| e.unit)
            .min()
            .unwrap_or(0)
    });
    debug!("Tetst cases are sorted.");
    let mut predicts = vec![];
    let mut chunks = construct_predictor(training, &alignments, contig)?;
    debug!("Predictors are constructed:{}", chunks.len());
    for (idx, (os, ms)) in chunks.iter().enumerate() {
        debug!("{}\t{}\t{}", idx, os.len(), ms.len());
    }
    let mut f = Factory::new();
    for read in tests {
        debug!(
            "{}:{}\t{}",
            read.id(),
            answer.contains_key(read.id()),
            dist.contains_key(read.id())
        );
        let id = read.id().to_string();
        let answer = if answer[&id] { 0 } else { 1 };
        let length = read.recover_raw_sequence().len();
        let dist = dist[&id];
        debug!("Read:{}\t{}", read, answer);
        let predict = make_prediction(&chunks, &read, &mut f);
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

type Seq = Vec<u8>;
// UnitID -> Sequences -> (color, seq)
fn construct_predictor(
    training: Vec<Record>,
    alignments: &[LastTAB],
    contig: &Contigs,
) -> Option<Vec<(Vec<Seq>, Vec<Seq>)>> {
    let is_original: HashMap<_, _> = training
        .iter()
        .map(|read| {
            let is_original = read
                .desc()
                .unwrap()
                .split(',')
                .nth(0)
                .unwrap()
                .contains("sample1");
            (read.id().to_string(), is_original)
        })
        .collect();
    let training = last_tiling::encoding(&training, contig, &alignments);
    debug!("Converted!");
    let mut chunks = vec![(vec![], vec![]); contig.get_last_unit(0)? as usize + 1];
    use last_tiling::revcmp;
    for read in training {
        for unit in read.seq().iter().filter_map(|e| e.encode()) {
            assert_eq!(0, unit.contig);
            let u = unit.unit as usize;
            let seq = if unit.is_forward {
                unit.bases.as_bytes().to_vec()
            } else {
                revcmp(unit.bases.as_bytes())
            };
            if is_original[read.id()] {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    Some(chunks)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    let (x1, x2) = ((x1 - log_denominator).exp(), (x2 - log_denominator).exp());
    assert!((x1 + x2 - 1.).abs() < 0.0001, "{}", x1 + x2);
    (x1, x2)
}

fn make_prediction(
    chunks: &[(Vec<Seq>, Vec<Seq>)],
    read: &last_tiling::EncodedRead,
    f: &mut Factory,
) -> u8 {
    let start = time::Instant::now();
    let predicts: Vec<(f64, f64)> = read
        .seq()
        .iter()
        .filter_map(|e| e.encode())
        .filter_map(|e| {
            let u = e.unit as usize;
            let (ref original, ref mutant) = &chunks[u];
            if original.len() < 5 || mutant.len() < 5 {
                return None;
            }
            let query = if e.is_forward {
                e.bases.as_bytes().to_vec()
            } else {
                last_tiling::revcmp(e.bases.as_bytes())
            };
            let o = unit_prediction(original, &query, f);
            let m = unit_prediction(mutant, &query, f);
            debug!(
                "Pred({})\t{}\t{}->{}\t{}",
                u,
                original.len(),
                mutant.len(),
                o,
                m
            );
            Some(as_weight(o, m))
        })
        .collect();
    if predicts.is_empty() {
        debug!("Overall:{:.4}\t{:.4}", 0.5, 0.5);
        return 0;
    }
    // Sum ln2 - H(p)
    for &(o, m) in &predicts {
        assert!(!o.is_nan() && !m.is_nan(), "{}\t{}", o, m);
    }
    let xlnx = |x: &f64| if x == &0. { 0. } else { x * x.ln() };
    let entropy: Vec<_> = predicts
        .iter()
        .map(|(o, m)| (2f64).ln() + xlnx(o) + xlnx(m))
        .collect();
    let sum = entropy.iter().sum::<f64>();
    assert!(!sum.is_nan() && sum > 0.0000001);
    let (p1, p2) = predicts
        .iter()
        .zip(entropy.iter())
        .map(|((o, m), w)| (o * w, m * w))
        .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
    let (p1, p2) = (p1 / sum, p2 / sum);
    assert!((p1 + p2 - 1.).abs() < 0.0001, "{}", p1 + p2);
    eprintln!("Overall:{:.4}\t{:.4}", p1, p2);
    debug!("Time {}", (time::Instant::now() - start).as_millis());
    if p2 < p1 {
        0
    } else {
        1
    }
}

fn unit_prediction(units: &[Seq], query: &[u8], f: &mut Factory) -> f64 {
    //let hmm = DBGHMM::new(units, K);
    let start = time::Instant::now();
    let hmm = f.generate(units, K);
    let res = hmm.forward(&query, &DEFAULT_CONFIG);
    assert!(!res.is_nan());
    debug!(
        "Unit {} {}",
        units.len(),
        (time::Instant::now() - start).as_millis()
    );
    res
}
