extern crate bio;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rayon;
extern crate kmer_counting;
use rayon::prelude::*;
fn main()->std::io::Result<()>{
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    debug!("Start process");
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = if args[1].ends_with('q'){
        bio::io::fastq::Reader::from_file(&args[1])?
        .records()
            .filter_map(|e| e.ok())
            .map(|e|e.seq().to_vec())
            .collect()
    }else{
        bio::io::fasta::Reader::from_file(&args[1])?
        .records()
            .filter_map(|e| e.ok())
            .map(|e|e.seq().to_vec())
            .collect()
    };
    let k: usize = args[2].parse().unwrap();
    debug!(
        "Opened records. {} bases.",
        reads.iter().map(|e| e.len()).sum::<usize>()
    );
    use std::io::{BufWriter, Write};
    let mut result:Vec<_> = kmer_counting::kmer_counting_native(&reads,k)
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
