#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate kmer_counting;
extern crate rayon;
use bio::io::fasta;
use kmer_counting::kmer_counting;
use rayon::prelude::*;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    debug!("Start process");
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = fasta::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let k: usize = args[2].parse().unwrap();
    debug!(
        "Opened fastq. {} bases.",
        reads.iter().map(|e| e.seq().len()).sum::<usize>()
    );
    use std::io::{BufWriter, Write};
    let reads: Vec<&[u8]> = reads.iter().map(|e| e.seq()).collect();
    let mut result: Vec<_> = kmer_counting(&reads, k)
        .into_iter()
        .filter(|(_, val)| val > &2)
        .collect();
    debug!("sorting");
    result.par_sort_by_key(|e| e.1);
    debug!("sorted");
    let out = std::io::stdout();
    let mut out = BufWriter::new(out.lock());
    debug!("start outputting");
    for (kmer, count) in result {
        let entropy = calc_entropy(&kmer.seq);
        out.write_all(&kmer.seq)?;
        out.write_all(b"\t")?;
        out.write_all(count.to_string().as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(entropy.to_string().as_bytes())?;
        out.write_all(b"\n")?;
    }
    debug!("End output");
    // use std::collections::HashMap;
    // let freqfreq:HashMap<_,i32> =
    //     result
    //         .into_iter()
    //         .fold(std::collections::HashMap::new(), |mut res, (_, val)| {
    //             let x = res.entry(val).or_default();
    //             *x += 1;
    //             res
    //         });
    // for (count, times) in freqfreq {
    //     writeln!(&mut out, "{}\t{}", count, times)?
    // }
    Ok(())
}

const K: usize = 3;
const MASK: usize = 0b111111;
const BUFSIZE: usize = 4 * 4 * 4;
fn calc_entropy(seq: &[u8]) -> f64 {
    let mut idx = seq[0..K - 1]
        .iter()
        .map(|&b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!(),
        })
        .fold(0, |x, y| x << 2 + y);
    let mut count: [u8; BUFSIZE] = [0; BUFSIZE];
    for b in &seq[K - 1..] {
        let next = match &b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!(),
        };
        idx = ((idx << 2) & MASK) + next;
        count[idx] += 1;
    }
    let tot = (seq.len() + 1 - K) as f64;
    -count.into_iter().fold(0., |x, &count| {
        if count == 0 {
            x
        } else {
            let f = count as f64 / tot;
            x + f * f.log2()
        }
    })
}
