extern crate bio_utils;
extern crate env_logger;
extern crate rayon;
#[macro_use]
extern crate log;
// use rayon::prelude::*;
const SD: f64 = 8.;
// use std::collections::HashMap;
use std::io::{BufWriter, Write};
// use std::sync::{Arc, Mutex};
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let mismatch_prob: f64 = args[2].parse().unwrap();
    let maf: Vec<_> = bio_utils::maf::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let subst: Vec<_> = maf.clone().into_iter().filter_map(Aln::new).collect();
    let ins: Vec<_> = maf.clone().into_iter().filter_map(Ins::new).collect();
    let del: Vec<_> = maf.clone().into_iter().filter_map(Del::new).collect();
    let mut pileups: Vec<_> = (0..1_000_000).map(PileUp::new).collect();
    let mut pileups_ins: Vec<_> = (0..1_000_000).map(PileUp::new).collect();
    let mut pileups_del: Vec<_> = (0..1_000_000).map(PileUp::new).collect();
    debug!("Converted {} alignments into minimal mode", subst.len());
    debug!("Register each alignment into pileups");
    for aln in &subst {
        for (idx, base) in aln.seq.iter().enumerate() {
            pileups
                .get_mut(idx + aln.start)
                .map(|e| e.push(*base))
                .unwrap();
        }
    }
    for aln in ins {
        for &(idx, base) in aln.seq.iter() {
            pileups_ins.get_mut(idx).map(|e| e.push(base)).unwrap();
        }
    }
    for aln in del {
        for &(idx, base) in aln.seq.iter() {
            pileups_del.get_mut(idx).map(|e| e.push(base)).unwrap();
        }
    }
    debug!("Write information on {}", &args[3]);
    let mut wtr = BufWriter::new(std::fs::File::create(&args[3])?);
    writeln!(&mut wtr, "Pos\tMAF\tCoverage\tType")?;
    for ((s, i), d) in pileups
        .iter()
        .zip(pileups_ins.iter())
        .zip(pileups_del.iter())
        .filter(|&((e, _), _)| e.coverage() > 20)
    {
        assert!(s.pos == i.pos && i.pos == d.pos);
        writeln!(&mut wtr, "{}\t{}\t{}\tSub", s.pos, s.maf(), s.coverage())?;
        let ins = i.coverage() as f64 / (s.coverage() + i.coverage()) as f64;
        writeln!(&mut wtr, "{}\t{}\t{}\tIns", i.pos, ins, i.coverage())?;
        let del = d.coverage() as f64 / (s.coverage() + d.coverage()) as f64;
        writeln!(&mut wtr, "{}\t{}\t{}\tDel", d.pos, del, d.coverage())?;
    }
    // debug!("Calling Variants");
    // let variants: Vec<_> = pileups
    //     .into_iter()
    //     .filter_map(|e| e.call(mismatch_prob))
    //     .collect();
    // debug!("{}", variants.len());
    // let stdout = std::io::stdout();
    // let mut wtr = BufWriter::new(stdout.lock());
    // let unwrap = |xs| match xs {
    //     Ok(res) => res,
    //     Err(why) => why,
    // };
    // debug!("Collecting spanning reads");
    // let result: Arc<Mutex<Vec<HashMap<_, (u16, u16, u16, u16, u16)>>>> =
    //     Arc::new(Mutex::new(vec![HashMap::new(); variants.len()]));
    // alns.par_iter().for_each(|aln| {
    //     let start = aln.start;
    //     let end = start + aln.seq.len();
    //     let start = unwrap(variants.binary_search_by_key(&start, |&(ref e, _)| e.pos));
    //     let end = unwrap(variants.binary_search_by_key(&end, |&(ref e, _)| e.pos));
    //     let mut temp = vec![];
    //     (start..end).into_iter().for_each(|i1| {
    //         let &(ref v1, base1) = &variants[i1];
    //         variants[i1 + 1..end].iter().for_each(|&(ref v2, base2)| {
    //             if aln.does_share(v1, v2) {
    //                 let is_minor1 = aln.is_minor(v1, base1);
    //                 let is_minor2 = aln.is_minor(v2, base2);
    //                 // I'm sure that v1.pos < v2.pos!
    //                 temp.push((i1, v2.pos, is_minor1, is_minor2));
    //             }
    //         });
    //     });
    //     let mut inner = result.lock().unwrap();
    //     for (i1, p2, is_minor1, is_minor2) in temp {
    //         let entry = inner.get_mut(i1).unwrap().entry(p2).or_default();
    //         entry.0 += 1;
    //         entry.1 += is_minor1 as u16;
    //         entry.2 += 1;
    //         entry.3 += is_minor2 as u16;
    //         entry.4 += (is_minor1 && is_minor2) as u16;
    //     }
    // });
    // debug!("Dump");
    // writeln!(&mut wtr, "pos1\tpos2\ttot1\tmac1\ttot2\tmac2\tshare")?;
    // let result: Vec<_> = Arc::try_unwrap(result).unwrap().into_inner().unwrap();
    // for (i1, rest) in result.into_iter().enumerate() {
    //     for (p2, (tot1, mac1, tot2, mac2, share)) in rest {
    //         if share > 3 {
    //             let p1 = variants[i1].0.pos;
    //             writeln!(
    //                 &mut wtr,
    //                 "{}\t{}\t{}\t{}\t{}\t{}\t{}",
    //                 p1, p2, tot1, mac1, tot2, mac2, share
    //             )?;
    //         }
    //     }
    // }
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
const BASES: &[u8; 5] = b"ACGT-";
impl PileUp {
    fn coverage(&self) -> usize {
        self.composition.iter().sum::<usize>()
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

// const BASES: [u8; 5] = *b"ACGT-";
impl Aln {
    fn new(aln: bio_utils::maf::Record) -> Option<Self> {
        if aln.score()? < 1_000. {
            return None;
        };
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
        Some(Self { seq, start })
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

#[derive(Debug, Clone)]
struct Ins {
    start: usize,
    seq: Vec<(usize, u8)>,
}

impl Ins {
    fn new(aln: bio_utils::maf::Record) -> Option<Self> {
        if aln.score()? < 1_000. {
            return None;
        };
        let start = aln.sequence()[0].start() as usize;
        let mut pos = start;
        let mut is_ins = false;
        let seq = aln.sequence()[0]
            .text()
            .iter()
            .zip(aln.sequence()[1].text().iter())
            .filter_map(|(&rfr, &qry)| match (rfr == b'-', is_ins) {
                (true, false) => {
                    is_ins = true;
                    Some((pos, qry.to_ascii_uppercase()))
                }
                (false, true) => {
                    is_ins = false;
                    pos += 1;
                    None
                }
                (true, true) => None,
                (false, false) => {
                    pos += 1;
                    None
                }
            })
            .collect();
        Some(Self { seq, start })
    }
}

#[derive(Debug, Clone)]
struct Del {
    start: usize,
    seq: Vec<(usize, u8)>,
}

impl Del {
    fn new(aln: bio_utils::maf::Record) -> Option<Self> {
        if aln.score()? < 1_000. {
            return None;
        };
        let start = aln.sequence()[0].start() as usize;
        let mut pos = start;
        let mut is_del = false;
        let seq = aln.sequence()[0]
            .text()
            .iter()
            .zip(aln.sequence()[1].text().iter())
            .filter_map(|(&rfr, &qry)| match (qry == b'-', is_del) {
                (true, false) => {
                    is_del = true;
                    pos += 1;
                    Some((pos - 1, rfr.to_ascii_uppercase()))
                }
                (false, true) => {
                    is_del = false;
                    pos += (rfr != b'-') as usize;
                    None
                }
                (true, true) => {
                    pos += 1;
                    None
                }
                (false, false) => {
                    is_del = false;
                    pos += (rfr != b'-') as usize;
                    None
                }
            })
            .collect();
        Some(Self { seq, start })
    }
}
