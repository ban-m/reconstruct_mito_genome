extern crate bio_utils;
extern crate rand;
use bio_utils::fasta;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use std::fs::File;
use std::io::BufWriter;
// 0.2%, 0.65%, 0.65%.
const SUB: f64 = 0.002;
const DEL: f64 = 0.0065;
const IN: f64 = 0.0065;

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let input = &fasta::parse_into_vec(&args[1])?[0];
    let seq: Vec<_> = input.seq().iter().map(|e| e.to_ascii_uppercase()).collect();
    let wtr = File::create("./data/mito_mutant.fasta")?;
    let mut wtr = fasta::Writer::new(BufWriter::new(wtr));
    let header = "depth=1.0 circular=true".to_string();
    wtr.write_record(&fasta::Record::with_data(input.id(), &Some(header), &seq))?;
    // seed.
    let seq = introduce_randomness(&seq, 43423);
    let header = "depth=0.1 circular=true".to_string();
    wtr.write_record(&fasta::Record::with_data("mutants", &Some(header), &seq))?;
    Ok(())
}

enum Op {
    Match,
    MisMatch,
    Del,
    In,
}

impl Op {
    fn weight(&self) -> f64 {
        match self {
            Op::Match => 1. - SUB - DEL - IN,
            Op::MisMatch => SUB,
            Op::Del => DEL,
            Op::In => IN,
        }
    }
}

const OPERATIONS: [Op; 4] = [Op::Match, Op::MisMatch, Op::Del, Op::In];

fn introduce_randomness(seq: &[u8], seed: u64) -> Vec<u8> {
    let mut res = vec![];
    let mut remainings: Vec<_> = seq.iter().copied().rev().collect();
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
    while !remainings.is_empty() {
        match OPERATIONS.choose_weighted(&mut rng, Op::weight).unwrap() {
            &Op::Match => res.push(remainings.pop().unwrap()),
            &Op::MisMatch => res.push(choose_base(&mut rng, remainings.pop().unwrap())),
            &Op::In => res.push(random_base(&mut rng)),
            &Op::Del => {
                remainings.pop().unwrap();
            }
        }
    }
    res
}

fn choose_base(rng: &mut StdRng, base: u8) -> u8 {
    let bases: Vec<u8> = b"ATCG".iter().filter(|&&e| e == base).copied().collect();
    *bases.choose_weighted(rng, |_| 1. / 3.).unwrap()
}

fn random_base(rng: &mut StdRng) -> u8 {
    *b"ATGC".choose_weighted(rng, |_| 1. / 4.).unwrap()
}
