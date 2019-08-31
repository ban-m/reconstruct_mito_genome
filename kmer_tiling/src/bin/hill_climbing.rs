extern crate bio;
extern crate handmade_bloom_filter;
extern crate rand;
use bio::io;
use handmade_bloom_filter::UpperBoundBFFactory;
use handmade_bloom_filter::HUGE_MODULO;
use std::io::{BufWriter, Write};
// use std::path::Path;
// use std::fs::File;
use rand::{seq::SliceRandom, thread_rng, Rng};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let t: usize = args[2].parse().unwrap();
    let input: Vec<_> = io::fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| e.seq().to_vec())
        .collect();
    eprintln!("Collect reads:{}", input.len());
    let bf = UpperBoundBFFactory::default()
        .k(20)
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    eprintln!("{}", bf);
    let mut rng = thread_rng();
    let mut wtr = BufWriter::new(std::io::stdout());
    //writeln!(&mut wtr, "ID\tIteration\tCount")?;
    for _t in 0..1 {
        let line = input.choose(&mut rng).unwrap();
        if line.len() >= 20 {
            let start = rng.gen_range(0, line.len() - 20);
            let kmer = line[start..start + 20].to_vec();
            let ub = bf.upper_bound(&kmer).unwrap();
            writeln!(&mut wtr, "{}\t{}", String::from_utf8_lossy(&kmer), ub)?;
            let candidate = candidate(&mut rng, 100, 20);
            for new_kmer in candidate {
                let new_ub = bf.upper_bound(&kmer).unwrap();
                let dist = edit_dist(&new_kmer, &kmer);
                writeln!(
                    &mut wtr,
                    "{}\t{}\t{}",
                    String::from_utf8_lossy(&new_kmer),
                    dist,
                    new_ub
                )?;
            }
        }
    }
    Ok(())
}

const NOTA: &[u8] = b"TCG";
const NOTC: &[u8] = b"ATG";
const NOTG: &[u8] = b"ATC";
const NOTT: &[u8] = b"ACG";
const BASE: &[u8] = b"ACGT";
#[allow(dead_code)]
fn mutate(kmer: &[u8], rng: &mut rand::rngs::ThreadRng) -> Vec<u8> {
    let mut kmer = kmer.to_vec();
    if let Some(point) = kmer.choose_mut(rng) {
        match *point { 
           b'A' => *point = *NOTA.choose(rng).unwrap(),
            b'C' => *point = *NOTC.choose(rng).unwrap(),
            b'G' => *point = *NOTG.choose(rng).unwrap(),
            b'T' => *point = *NOTT.choose(rng).unwrap(),
            _ => panic!(),
        }
    }
    kmer
}

type Kmer = Vec<u8>;
// generate t random kmer
fn candidate(rng: &mut rand::rngs::ThreadRng, t: usize, k: usize) -> Vec<Kmer> {
    (0..t).map(|_| generate(rng, k)).collect()
}

fn generate(rng: &mut rand::rngs::ThreadRng, k: usize) -> Kmer {
    (0..k)
        .filter_map(|_| BASE.choose(rng))
        .map(|e| *e)
        .collect()
}

fn edit_dist(x: &[u8], y: &[u8]) -> usize {
    let mut d = vec![vec![0; y.len() + 1]; x.len() + 1];
    for i in 0..=x.len() {
        d[i][0] = i;
    }
    for j in 0..=y.len() {
        d[0][j] = j;
    }
    for i in 1..=x.len() {
        for j in 1..=y.len() {
            let is_match = if x[i-1] == y[j-1] { 0 } else { 1 };
            d[i][j] = (d[i - 1][j - 1] + is_match)
                .min(d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1);
        }
    }
    d[x.len()][y.len()]
}
