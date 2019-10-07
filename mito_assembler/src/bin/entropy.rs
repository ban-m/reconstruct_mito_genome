// calc k-mer entropy on each read.
extern crate bio_utils;
extern crate rayon;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::ops::Shl;
fn main() -> std::io::Result<()> {
    let reads = std::env::args()
        .nth(1)
        .map(|e| bio_utils::fasta::parse_into_vec(e))
        .unwrap()?;
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    // stdout.write_all(b"k\tentropy\n")?;
    // let ks = vec![3, 4, 5, 6];
    // for read in reads {
    //     for k in &ks {
    //         let entropy = calc_entropy(read.seq(), *k);
    //         writeln!(&mut stdout, "{}\t{}", k, entropy)?;
    //     }
    // }
    let reads: Vec<_> = reads
        .into_par_iter()
        .filter(|read| {
            let e = calc_entropy(read.seq(), 6);
            8.0 < e && e < 9.0
        })
        .collect();
    for read in reads {
        writeln!(&mut stdout, "{}", read)?;
    }
    stdout.flush()
}

fn calc_entropy(read: &[u8], k: usize) -> f64 {
    if read.len() < k {
        0.
    } else {
        let mut slots: Vec<u32> = vec![0; 4usize.pow(k as u32)];
        let mask = 4usize.pow(k as u32) - 1;
        let mut current = calc_index(&read[..k - 1]);
        let total = (read.len() - k + 1) as f64;
        for base in &read[k..] {
            current = current.shl(2) & mask;
            current += match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => unreachable!(),
            };
            slots[current] += 1;
        }
        slots
            .into_iter()
            .filter(|&count| count != 0)
            .map(|e| e as f64 / total)
            .map(|e| -e * e.log2())
            .sum::<f64>()
    }
}

fn calc_index(seq: &[u8]) -> usize {
    seq.iter()
        .map(|base| match base {
            b'A' | b'a' => 0usize,
            b'C' | b'c' => 1usize,
            b'G' | b'g' => 2usize,
            b'T' | b't' => 3usize,
            _ => unreachable!(),
        })
        .fold(0, |sum, b| sum.shl(2) + b)
}
