extern crate bio_utils;
extern crate last_tiling;
extern crate rayon;
use rayon::prelude::*;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let tab: Vec<_> = last_tiling::parse_tab_file(&args[1])?;
    let read: Vec<_> = bio_utils::fasta::parse_into_vec(&args[2])?;
    let contig: Vec<_> = bio_utils::fasta::parse_into_vec(&args[3])?;
    eprintln!("{}", tab.len());
    println!("Pm\tPi\tPd");
    let summary = summarize_tab(tab, read, contig);
    println!("{}", summary);
    let mismatch_prob = summary.p_mismatch;
    println!("{},{}", mismatch_prob, mismatch_prob / 3.);
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct Summary {
    p_mismatch: f64,
    p_in: f64,
    p_del: f64,
}

impl Summary {
    fn new(p_mismatch: f64, p_in: f64, p_del: f64) -> Self {
        Self {
            p_mismatch,
            p_in,
            p_del,
        }
    }
}

impl std::fmt::Display for Summary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:.2}\t{:.2}\t{:.2}",
            100. * self.p_mismatch,
            100. * self.p_in,
            100. * self.p_del,
        )
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Op {
    Match,
    Mism,
    Del,
    In,
}

fn summarize_tab(
    tab: Vec<last_tiling::LastTAB>,
    read: Vec<bio_utils::fasta::Record>,
    contig: Vec<bio_utils::fasta::Record>,
) -> Summary {
    let reads: HashMap<_, _> = read.into_iter().map(|e| (e.id().to_string(), e)).collect();
    let contigs: HashMap<_, _> = contig
        .into_iter()
        .map(|e| (e.id().to_string(), e))
        .collect();
    let ops: Vec<Vec<_>> = tab
        .into_par_iter()
        .map(|tab| {
            let rn = tab.seq2_name();
            let cn = tab.seq1_name();
            recover(&reads[rn], &contigs[cn], tab)
        })
        .collect();
    summarize_operations(ops)
}

fn recover(
    query: &bio_utils::fasta::Record,
    refr: &bio_utils::fasta::Record,
    tab: last_tiling::LastTAB,
) -> Vec<Op> {
    let (mut r, mut q) = (tab.seq1_start(), tab.seq2_start());
    let refr = refr.seq();
    let query = if tab.seq2_direction().is_forward() {
        query.seq().to_vec()
    } else {
        last_tiling::revcmp(query.seq())
    };
    let mut ops = vec![];
    for op in tab.alignment() {
        use last_tiling::Op::*;
        match op {
            Match(l) => {
                ops.extend(
                    refr[r..r + l]
                        .iter()
                        .zip(query[q..q + l].iter())
                        .map(|(a, b)| if a == b { Op::Match } else { Op::Mism }),
                );
                r += l;
                q += l;
            }
            Seq1In(l) => {
                ops.extend(vec![Op::In; l]);
                q += l;
            }
            Seq2In(l) => {
                ops.extend(vec![Op::Del; l]);
                r += l;
            }
        }
    }
    ops
}

fn summarize_operations(opss: Vec<Vec<Op>>) -> Summary {
    let mut mismatch = 0;
    let mut num_del = 0;
    let mut num_in = 0;
    use Op::*;
    let total:usize = opss.iter().map(|e| e.len()).sum();
    for ops in opss {
        for op in ops {
            match op {
                Mism => mismatch += 1,
                Del => num_del += 1,
                In => num_in += 1,
                _ => {},
            }
        }
    }
    let div = |x, y| x as f64 / y as f64;
    //eprintln!("{}\t{}\t{}", matchmis, num_in,num_del);
    Summary::new(
        div(mismatch, total),
        div(num_in, total),
        div(num_del, total),
    )
}
