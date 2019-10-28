extern crate bio_utils;
extern crate env_logger;
extern crate rayon;
#[macro_use]
extern crate log;
use rayon::prelude::*;
const SD: f64 = 8.;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let maf: Vec<_> = bio_utils::maf::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let mismatch_prob: f64 = args[2].parse().unwrap();
    debug!("Convert alignment into minimal mode");
    let alns: Vec<_> = maf.into_iter().map(Aln::new).collect();
    let mut pileups: Vec<_> = (0..1_000_000).map(PileUp::new).collect();
    debug!("Register each alignment into pileups");
    for aln in &alns {
        for (idx, base) in aln.seq.iter().enumerate() {
            pileups
                .get_mut(idx + aln.start)
                .map(|e| e.push(*base))
                .unwrap();
        }
    }
    debug!("Write information on ./result/maf.tsv");
    {
        let mut wtr = BufWriter::new(std::fs::File::create("./result/maf.tsv")?);
        for v in pileups.iter().filter(|&e| e.nonzero()) {
            writeln!(&mut wtr, "{}\t{}", v.pos, v.maf())?;
        }
    }
    debug!("Calling Variants");
    let variants: Vec<_> = pileups
        .into_iter()
        .filter_map(|e| e.call(mismatch_prob))
        .collect();
    debug!("{}", variants.len());
    use std::io::{BufWriter, Write};
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    let result = {
        debug!("Collecting spanning reads");
        let counts: Vec<_> = alns
            .par_iter()
            .map(|aln| {
                let mut result: HashMap<_, (i32, i32, i32, i32, i32)> = HashMap::new();
                let start = aln.start;
                let end = start + aln.seq.len();
                let start = match variants.binary_search_by_key(&start, |&(ref e, _)| e.pos) {
                    Ok(res) => res,
                    Err(why) => why,
                };
                let end = match variants.binary_search_by_key(&end, |&(ref e, _)| e.pos) {
                    Ok(res) => res,
                    Err(why) => why,
                };
                (start..end).into_iter().for_each(|i1| {
                    let &(ref v1, base1) = &variants[i1];
                    (i1 + 1..end).for_each(|i2| {
                        // Count the link strength between i-th and j-th variants.
                        let &(ref v2, base2) = &variants[i2];
                        if aln.does_share(v1, v2) {
                            let is_minor1 = aln.is_minor(v1, base1);
                            let is_minor2 = aln.is_minor(v2, base2);
                            let entry = result.entry((v1.pos, v2.pos)).or_default();
                            entry.0 += 1;
                            entry.1 += is_minor1 as i32;
                            entry.2 += 1;
                            entry.3 += is_minor2 as i32;
                            entry.4 += (is_minor1 && is_minor2) as i32;
                        }
                    })
                });
                result
            })
            .collect();
        let mut result = HashMap::new();
        for block in counts {
            for (key, (x0, x1, x2, x3, x4)) in block {
                let entry = result.entry(key).or_insert((0, 0, 0, 0, 0));
                entry.0 += x0;
                entry.1 += x1;
                entry.2 += x2;
                entry.3 += x3;
                entry.4 += x4;
            }
        }
        result
    };
    debug!("Dump");
    writeln!(&mut wtr, "pos1\tpos2\tmac1\ttot1\tmac2\ttot2\tshare")?;
    for ((p1, p2), (mac1, tot1, mac2, tot2, share)) in result {
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p1, p2, mac1, tot1, mac2, tot2, share
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct PileUp {
    pos: usize,
    // A,C,G,T,gap.
    composition: [usize; 5],
}

impl std::fmt::Display for PileUp {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Pos:{}\tA:{}\tC:{}\tG:{}\tT:{}\tGap:{}",
            self.pos,
            self.composition[0],
            self.composition[1],
            self.composition[2],
            self.composition[3],
            self.composition[4]
        )
    }
}

impl PileUp {
    fn nonzero(&self) -> bool {
        self.composition.iter().sum::<usize>() != 0
    }
    fn push(&mut self, base: u8) {
        match base {
            b'A' => self.composition[0] += 1,
            b'C' => self.composition[1] += 1,
            b'G' => self.composition[2] += 1,
            b'T' => self.composition[3] += 1,
            b'-' => {} //Disabled currently  self.composition[4] += 1,
            _ => unreachable!(),
        }
    }
    fn maf(&self) -> f64 {
        let (idx, _) = self
            .composition
            .iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();
        let (_, second_major) = self
            .composition
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .max_by_key(|(_, &c)| c)
            .unwrap();
        let total = self.composition.iter().sum::<usize>() as f64;
        let second_major = *second_major as f64;
        second_major / total
    }
    fn call(self, p: f64) -> Option<(Self, u8)> {
        let (idx, _) = self
            .composition
            .iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();
        let (second_major_idx, second_major) = self
            .composition
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .max_by_key(|(_, &c)| c)
            .unwrap();
        let base = BASES[second_major_idx];
        let total = self.composition.iter().sum::<usize>() as f64;
        let second_major = *second_major as f64;
        let sd = (p * (1. - p) / total).sqrt();
        let thr = p + SD * sd;
        if second_major / total > thr {
            debug!("{}\t{}\t{}", self, thr, second_major / total);
            Some((self, base))
        } else {
            None
        }
    }
    fn new(pos: usize) -> Self {
        let composition = [0, 0, 0, 0, 0];
        Self { pos, composition }
    }
}

#[derive(Debug, Clone)]
struct Aln {
    start: usize,
    seq: Vec<u8>,
}

const BASES: [u8; 5] = *b"ACGT-";
impl Aln {
    fn new(aln: bio_utils::maf::Record) -> Self {
        let start = aln.sequence()[0].start() as usize;
        let non_gap = b"ACGT";
        let seq = aln.sequence()[0]
            .text()
            .iter()
            .zip(aln.sequence()[1].text().iter())
            .filter_map(|(&rfr, &qry)| {
                if non_gap.contains(&rfr) {
                    Some(qry)
                } else {
                    None
                }
            })
            .map(|e| e.to_ascii_uppercase())
            .collect();
        Self { seq, start }
    }
    fn does_share(&self, v1: &PileUp, v2: &PileUp) -> bool {
        // Include check.
        let end = self.start + self.seq.len();
        self.start <= v1.pos && v1.pos < end && self.start <= v2.pos && v2.pos < end
    }
    fn is_minor(&self, v: &PileUp, b: u8) -> bool {
        self.seq[v.pos - self.start] == b
    }
}
