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
const K: usize = 6;
// const MAX_CHAIN: usize = 10;
use dbg_hmm::*;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("predict_mockdata=debug"))
        .init();
    info!("PRED\tReadID\tPrediction\tIsCorrect\tLKDiff");
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
    let (training, testset): (Vec<_>, Vec<_>) = reads.into_iter().partition(|_| rng.gen_bool(0.1));
    // let (training, testset): (Vec<_>, Vec<_>) =
    //     reads.into_iter().partition(|r| dist[r.id()] < len / 2);
    debug!(
        "{} training reads and {} test reads dataset.",
        training.len(),
        testset.len()
    );
    let result = predict(training, testset, alignments, &reference, &answer).unwrap();
    debug!("Dump results");
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    writeln!(
        &mut wtr,
        "ReadID\tAnswer\tPredict\tDistance\tLength\tLKDiff"
    )?;
    let (mut t_p, mut t_n, mut f_p, mut f_n) = (0, 0, 0, 0);
    for (readid, predict, length, lk) in result {
        let answer = if answer[&readid] { 0 } else { 1 };
        let dist = dist[&readid];
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}\t{}",
            readid, answer, predict, dist, length, lk,
        )?;
        match (predict, answer) {
            (0, 0) => t_p += 1,
            (0, 1) => f_n += 1,
            (1, 0) => f_p += 1,
            (1, 1) => t_n += 1,
            _ => {}
        }
    }
    let sum = t_p + t_n + f_n + f_p;
    let ok = t_p + t_n;
    let acc = ok as f64 / sum as f64;
    let sens = t_p as f64 / (t_p + f_n) as f64;
    let spec = t_p as f64 / (t_p + f_p) as f64;
    info!("SUMMARY\tResult:{}\t{}\t{:.3}", sum, ok, acc);
    info!("SUMMARY\tResult:{:.3}\t{:.3}", sens, spec);
    Ok(())
}

fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
    ans: &HashMap<String, bool>,
) -> Option<Vec<(String, u8, usize, f64)>> {
    // let is_test: std::collections::HashSet<_> = tests.iter().map(|e| e.id().to_string()).collect();
    let (is_original, data, border) = setup(training, tests, alignments, contig);
    let len = contig.get_last_unit(0)? as usize + 1;
    let answer: Vec<_> = data.iter().map(|e| ans[e.id()]).collect::<Vec<_>>();
    // let lk = compute_likelihood(&data, &answer, len)
    //     .into_iter()
    //     .skip(border)
    //     .sum::<f64>();
    // info!("SUMMARY\tObjLK\tlog lk\t{:.2}", lk);
    let pred = initial_pred(&is_original, &data, border, len);
    // assert_eq!(lks.len(), data.len());
    // let lks = compute_likelihood_diff(&data, &pred, len);
    assert_eq!(pred.len(), data.len());
    let acc = pred
        .iter()
        .zip(answer.iter())
        .skip(border)
        .filter(|(f, g)| f == g)
        .count() as f64;
    info!(
        "Initial prediction's accuracy:{:.4}",
        acc / (pred.len() - border) as f64
    );
    let em_pred = em_prediction(&pred, &data, border, len);
    let lks = compute_likelihood_diff(&data, &em_pred, len);
    let result: Vec<_> = data
        .iter()
        .zip(em_pred.iter())
        // .zip(pred.iter())
        .zip(lks.into_iter())
        .skip(border)
        .map(|((read, predict), lk)| {
            let id = read.id().to_string();
            let pred = if *predict { 0u8 } else { 1u8 };
            let len = read.seq().iter().map(|e| e.bases().len()).sum::<usize>();
            (id, pred, len, lk)
        })
        .collect();
    Some(result)
}

fn setup(
    training: Vec<Record>,
    tests: Vec<Record>,
    alns: Vec<LastTAB>,
    contig: &Contigs,
) -> (Vec<bool>, Vec<ERead>, usize) {
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
    let training_size = training.len();
    training.append(&mut tests);
    let data: Vec<_> = training.into_iter().map(ERead::new).collect();
    (is_original, data, training_size)
}

#[allow(dead_code)]
fn compute_likelihood(data: &[ERead], is_original: &[bool], len: usize) -> Vec<f64> {
    let chunks = construct_predictors_simple(data, is_original, len);
    is_original
        .par_iter()
        .zip(data.par_iter())
        .map(|(is_ori, read)| {
            read.seq()
                .iter()
                .filter_map(|e| {
                    let u = e.unit as usize;
                    let ref refs = if *is_ori { &chunks[u].0 } else { &chunks[u].1 };
                    if refs.len() < 10 {
                        return None;
                    }
                    let query = e.bases();
                    Some(unit_predict(&query, &refs, K))
                })
                .sum::<f64>()
        })
        .collect()
}

fn compute_likelihood_diff(data: &[ERead], is_original: &[bool], len: usize) -> Vec<f64> {
    let chunks = construct_predictors_simple(data, is_original, len);
    is_original
        .par_iter()
        .zip(data.par_iter())
        .map(|(is_ori, read)| {
            read.seq()
                .iter()
                .filter_map(|e| {
                    let u = e.unit as usize;
                    let (ref ori, ref muta) = &chunks[u];
                    if ori.len() < 10 || muta.len() < 10 {
                        return None;
                    }
                    let query = e.bases();
                    let o = unit_predict(&query, &ori, K);
                    let m = unit_predict(&query, &muta, K);
                    let diff = if *is_ori { o - m } else { m - o };
                    Some(diff)
                })
                .sum::<f64>()
        })
        .collect()
}

fn initial_pred(init: &[bool], data: &[ERead], border: usize, len: usize) -> Vec<bool> {
    let chunks = construct_predictors_simple(&data[..border], &init[..border], len);
    let mut pred: Vec<_> = init[..border]
        .par_iter()
        .copied()
        .chain(data.par_iter().skip(border).map(|read| {
            let (predict, _lk) = make_prediction_simple(&chunks, &read);
            predict
        }))
        .collect();
    assert_eq!(data.len(), pred.len());
    let mut updated = true;
    while updated {
        let chunks = construct_predictors_simple(&data, &pred, len);
        updated = pred
            .par_iter_mut()
            .enumerate()
            .map(|(idx, p)| {
                let (predict, _lk) = make_prediction_simple(&chunks, &data[idx]);
                let up = *p != predict;
                *p = predict;
                up
            })
            .reduce(|| false, |p, q| p | q);
    }
    pred
}

fn em_prediction(init: &[bool], data: &[ERead], border: usize, len: usize) -> Vec<bool> {
    let (mut w0, mut w1) = {
        let len = init.len() as f64;
        let w0 = init.iter().filter(|&&e| e).count() as f64 / len;
        let w1 = init.iter().filter(|&&e| !e).count() as f64 / len;
        (w0, w1)
    };
    let mut gamma0: Vec<_> = init.iter().map(|&e| if e { 1. } else { 0. }).collect();
    let mut gamma1: Vec<_> = init.iter().map(|&e| if !e { 1. } else { 0. }).collect();
    let mut model0: Vec<_> = construct_predictors(data, &gamma0, len, K);
    let mut model1: Vec<_> = construct_predictors(data, &gamma1, len, K);
    let mut lk = compute_lk_from(&data[border..], &model0, &model1, w0, w1);
    info!("LK\t{:.4}\t0", lk);
    for i in 1..=10 {
        data.par_iter()
            .zip(gamma0.par_iter_mut())
            .zip(gamma1.par_iter_mut())
            .skip(border)
            .for_each(|((read, g0), g1)| {
                let log_m0 = read
                    .seq()
                    .iter()
                    .map(|s| model0[s.unit as usize].forward(s.bases(), &DEFAULT_CONFIG))
                    .sum::<f64>();
                let log_m1 = read
                    .seq()
                    .iter()
                    .map(|s| model1[s.unit as usize].forward(s.bases(), &DEFAULT_CONFIG))
                    .sum::<f64>();
                let w = calc_logsum(log_m0, log_m1, w0, w1);
                *g0 = (w0.ln() + log_m0 - w).exp();
                *g1 = (w1.ln() + log_m1 - w).exp();
            });
        let tot = gamma0.iter().sum::<f64>() + gamma1.iter().sum::<f64>();
        w0 = gamma0.iter().sum::<f64>() / tot;
        w1 = gamma1.iter().sum::<f64>() / tot;
        model0 = construct_predictors(&data, &gamma0, len, K);
        model1 = construct_predictors(&data, &gamma1, len, K);
        let next_lk = compute_lk_from(&data[border..], &model0, &model1, w0, w1);
        info!("LK\t{:.4}\t{}", lk, i);
        // debug!("(w0,w1)=({},{})", w0, w1);
        if (next_lk - lk).abs() < 0.001 {
            break;
        } else {
            lk = next_lk;
        }
    }
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1893744392);
    gamma0.iter().map(|&w| rng.gen_range(0., 1.) < w).collect()
}

fn make_prediction_simple(chunks: &[(Vec<&[u8]>, Vec<&[u8]>)], read: &ERead) -> (bool, f64) {
    let predicts: Vec<(f64, f64)> = read
        .seq()
        .iter()
        .filter_map(|e| {
            let u = e.unit as usize;
            let (ref original, ref mutant) = &chunks[u];
            if original.is_empty() || mutant.is_empty() {
                return None;
            }
            let query = e.bases();
            let min = original.len().min(mutant.len());
            let o = unit_predict(&query, &original[..min], K);
            let m = unit_predict(&query, &mutant[..min], K);
            Some((o, m))
        })
        // .take(MAX_CHAIN)
        .collect();
    let p1_minus_p2 = predicts.iter().map(|(l1, l2)| l1 - l2).sum::<f64>();
    if predicts.is_empty() {
        return (true, 0.);
    }
    (p1_minus_p2.is_sign_positive(), p1_minus_p2)
}

fn construct_predictors_simple<'a>(
    data: &'a [ERead],
    is_original: &[bool],
    len: usize,
) -> Vec<(Vec<&'a [u8]>, Vec<&'a [u8]>)> {
    let mut chunks = vec![(vec![], vec![]); len];
    for (&is_original, read) in is_original.iter().zip(data.iter()) {
        for unit in read.seq().iter() {
            let u = unit.unit as usize;
            let seq = unit.bases();
            if is_original {
                chunks[u].0.push(seq);
            } else {
                chunks[u].1.push(seq);
            }
        }
    }
    chunks
}

fn construct_predictors(data: &[ERead], weights: &[f64], len: usize, k: usize) -> Vec<DBGHMM> {
    let mut buf: Vec<Vec<&[u8]>> = vec![vec![]; len];
    let mut ws: Vec<Vec<f64>> = vec![vec![]; len];
    for (read, &w) in data
        .iter()
        .zip(weights.iter())
        .filter(|(_, &w)| w > 0.00001)
    {
        for c in read.seq().iter() {
            buf[c.unit as usize].push(c.bases());
            ws[c.unit as usize].push(w);
        }
    }
    let mut f = dbg_hmm::Factory::new();
    buf.into_iter()
        .zip(ws.into_iter())
        .map(|(cs, ws)| f.generate_with_weight(&cs, &ws, k))
        .collect()
}

fn compute_lk_from(data: &[ERead], m0: &[DBGHMM], m1: &[DBGHMM], w0: f64, w1: f64) -> f64 {
    data.par_iter()
        .map(|read| {
            read.seq()
                .iter()
                .map(|c| {
                    let m0 = m0[c.unit as usize].forward(c.bases(), &DEFAULT_CONFIG);
                    let m1 = m1[c.unit as usize].forward(c.bases(), &DEFAULT_CONFIG);
                    calc_logsum(m0, m1, w0, w1)
                })
                .sum::<f64>()
        })
        .sum::<f64>()
}
