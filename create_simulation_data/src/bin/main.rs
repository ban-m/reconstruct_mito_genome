extern crate bio_utils;
extern crate rand;
use bio_utils::fasta;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use std::fs::File;
use std::io::BufWriter;
// Repeat one
const REP1: (usize, usize) = (4193, 199_390);
const REP2: (usize, usize) = (56_409, 258_637);
// 0.1%, 0.45%, 0.45%.
const SUB: f64 = 0.003;
const DEL: f64 = 0.003;
const IN: f64 = 0.003;

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let input = &fasta::parse_into_vec(&args[1])?[0];
    std::fs::create_dir_all("./data/flip_repeats")?;
    let seq: Vec<_> = input.seq().iter().map(|e| e.to_ascii_uppercase()).collect();
    {
        let wtr = File::create("./data/flip_repeats/forward_repeat.fasta")?;
        let mut wtr = fasta::Writer::new(BufWriter::new(wtr));
        let header = "depth=1.0 circular=true".to_string();
        wtr.write_record(&fasta::Record::with_data("original", &Some(header), &seq))?;
        let (c1, c2) = split_forward(&seq);
        // seed.
        let c1 = introduce_randomness(&c1, 3424);
        let c2 = introduce_randomness(&c2, 43423);
        let header = "depth=0.5 circular=true".to_string();
        wtr.write_record(&fasta::Record::with_data(
            "right_child",
            &Some(header.clone()),
            &c1,
        ))?;
        wtr.write_record(&fasta::Record::with_data(
            "left_child",
            &Some(header.clone()),
            &c2,
        ))?;
    }
    {
        let wtr = File::create("./data/flip_repeats/reverse_repeat.fasta")?;
        let mut wtr = fasta::Writer::new(BufWriter::new(wtr));
        let header = "depth=1.0 circular=true".to_string();
        wtr.write_record(&fasta::Record::with_data("original", &Some(header), &seq))?;
        let c1 = split_reverse(&seq);
        // seed.
        let c1 = introduce_randomness(&c1, 3423424);
        let header = "depth=0.4 circular=true".to_string();
        wtr.write_record(&fasta::Record::with_data(
            "right_child",
            &Some(header.clone()),
            &c1,
        ))?;
    }
    Ok(())
}

fn split_forward(seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let c1: Vec<_> = seq[..REP1.0]
        .iter()
        .chain(seq[REP1.1..].iter())
        .copied()
        .collect();
    let c2: Vec<_> = seq[REP1.0..REP1.1].to_vec();
    (c1, c2)
}

fn split_reverse(seq: &[u8]) -> Vec<u8> {
    seq[..REP2.0]
        .iter()
        .copied()
        .chain(seq[REP2.0..REP2.1].iter().rev().map(|e| comp(*e)))
        .chain(seq[REP2.1..].iter().copied())
        .collect()
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => unreachable!(),
    }
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
