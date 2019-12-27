extern crate bio_utils;
extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta::Record;
use last_decompose::clustering;
use last_tiling::Contigs;
use last_tiling::LastTAB;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
const K: usize = 4;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = bio_utils::fasta::parse_into_vec(&args[1])?
        .into_iter()
        .filter(|r| r.desc().unwrap().contains("sample"))
        .collect();
    debug!("{} reads in total.", reads.len());
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[3])?;
    let config = {
        let contigs = bio_utils::fasta::parse_into_vec(&args[2]).unwrap();
        last_decompose::error_profile::summarize_tab(&alignments, &reads, &contigs)
    };
    debug!("{:?}", config);
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
    let (training, testset): (Vec<_>, Vec<_>) = if args[4] == "random" {
        debug!("Random mode");
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1893749823);
        let (training, testset): (Vec<_>, Vec<_>) =
            reads.into_iter().partition(|_| rng.gen_bool(0.27));
        (training, testset)
    } else {
        debug!("Half mode");
        reads.into_iter().partition(|r| dist[r.id()] < len / 2)
    };
    debug!(
        "{} training reads and {} test reads dataset.",
        training.len(),
        testset.len()
    );
    let s = std::time::Instant::now();
    let result = predict(training, testset, alignments, &reference, &answer, &config).unwrap();
    debug!("Elapsed Time:{:?}", std::time::Instant::now() - s);
    debug!("Dump results");
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    writeln!(&mut wtr, "ReadID\tAnswer\tPredict\tDistance")?;
    let (mut t_p, mut t_n, mut f_p, mut f_n) = (0, 0, 0, 0);
    for (readid, predict) in result {
        let answer = if answer[&readid] { 0 } else { 1 };
        let dist = dist[&readid];
        writeln!(&mut wtr, "{}\t{}\t{}\t{}", readid, answer, predict, dist,)?;
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
    info!("SUMMARY\tResult\tSens:Spec={:.3}:{:.3}", sens, spec);
    Ok(())
}

fn predict(
    training: Vec<Record>,
    tests: Vec<Record>,
    alignments: Vec<LastTAB>,
    contig: &Contigs,
    ans: &HashMap<String, bool>,
    config: &dbg_hmm::Config,
) -> Option<Vec<(String, u8)>> {
    let (data, border) = setup(training, tests, alignments, contig);
    let answer: Vec<_> = data.iter().map(|e| ans[e.id()]).collect::<Vec<_>>();
    let label: Vec<u8> = answer[..border]
        .iter()
        .map(|&e| if e { 0 } else { 1 })
        .collect();
    let forbid: Vec<_> = data.iter().map(|_| vec![]).collect();
    let contigs: Vec<_> = (0..contig.get_num_of_contigs())
        .map(|e| contig.get_last_unit(e as u16).unwrap() as usize + 1)
        .collect();
    debug!("{:?}", contigs);
    let result = {
        let answer: Vec<_> = answer
            .iter()
            .skip(border)
            .map(|&e| if e { 0 } else { 1 })
            .collect();
        clustering(&data, &label, &forbid, K, 2, &contigs, &answer, config)
    };
    let result: Vec<_> = data[border..]
        .iter()
        .zip(result)
        .map(|(read, cluster)| {
            let id = read.id().to_string();
            (id, cluster)
        })
        .collect();
    Some(result)
}

fn setup(
    training: Vec<Record>,
    tests: Vec<Record>,
    alns: Vec<LastTAB>,
    contig: &Contigs,
) -> (Vec<last_decompose::ERead>, usize) {
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
    let data: Vec<_> = training
        .into_iter()
        .map(last_decompose::ERead::new_no_gapfill)
        .collect();
    (data, training_size)
}
