extern crate bio_utils;
extern crate last_tiling;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let maf: Vec<_> = bio_utils::maf::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let tab: Vec<_> = last_tiling::parse_tab_file(&args[2])?;
    let read: Vec<_> = bio_utils::fasta::parse_into_vec(&args[3])?;
    let contig: Vec<_> = bio_utils::fasta::parse_into_vec(&args[4])?;
    eprintln!("{}", maf.len());
    eprintln!("{}", tab.len());
    // Sanitization check.
    println!("Pm\tPSi\tPSd\tPEi\tPEd\tPid\tPdi");
    {
        let summary = summarize_tab(tab, read, contig);
        println!("{}", summary);
    }
    {
        let summary = summarize_maf(maf);
        println!("{}", summary);
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct Summary {
    p_match: f64,
    p_start_in: f64,
    p_start_del: f64,
    p_ext_in: f64,
    p_ext_del: f64,
    p_in_to_del: f64,
    p_del_to_in: f64,
}

impl Summary {
    fn new(
        p_match: f64,
        p_start_in: f64,
        p_start_del: f64,
        p_ext_in: f64,
        p_ext_del: f64,
        p_in_to_del: f64,
        p_del_to_in: f64,
    ) -> Self {
        Self {
            p_match,
            p_start_in,
            p_start_del,
            p_ext_in,
            p_ext_del,
            p_in_to_del,
            p_del_to_in,
        }
    }
}

impl std::fmt::Display for Summary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
            100. * self.p_match,
            100. * self.p_start_in,
            100. * self.p_start_del,
            100. * self.p_ext_in,
            100. * self.p_ext_del,
            100. * self.p_in_to_del,
            100. * self.p_del_to_in
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
        .into_iter()
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

fn summarize_maf(maf: Vec<bio_utils::maf::Record>) -> Summary {
    let ops: Vec<Vec<_>> = maf.into_iter().map(maf_to_ops).collect();
    //eprintln!("{}",ops.iter().map(|e|e.len()).sum::<usize>());
    summarize_operations(ops)
}

fn maf_to_ops(maf: bio_utils::maf::Record) -> Vec<Op> {
    let seq = maf.sequence();
    seq[0]
        .text()
        .iter()
        .zip(seq[1].text().iter())
        .map(|(a, b)| {
            if a == b {
                Op::Match
            } else if a == &b'-' {
                Op::In
            } else if b == &b'-' {
                Op::Del
            } else {
                Op::Mism
            }
        })
        .collect()
}

fn summarize_operations(opss: Vec<Vec<Op>>) -> Summary {
    let mut matchmis = 0;
    let mut num_seq = 0;
    let mut num_del = 0;
    let mut num_in = 0;
    let mut mm_after_mm = 0;
    let mut in_after_mm = 0;
    let mut in_after_del = 0;
    let mut in_after_in = 0;
    let mut del_after_mm = 0;
    let mut del_after_del = 0;
    let mut del_after_in = 0;
    use Op::*;
    for ops in opss {
        num_seq += 1;
        matchmis += ops
            .iter()
            .filter(|&e| match e {
                Match | Mism => true,
                _ => false,
            })
            .count();
        num_del += ops
            .iter()
            .filter(|&e| match e {
                Del => true,
                _ => false,
            })
            .count();
        num_in += ops
            .iter()
            .filter(|&e| match e {
                In => true,
                _ => false,
            })
            .count();
        for before_after in ops.windows(2) {
            let b = before_after[0];
            let a = before_after[1];
            match (b, a) {
                (Match, Mism) | (Match, Match) | (Mism, Match) | (Mism, Mism) => mm_after_mm += 1,
                (Mism, Del) | (Match, Del) => del_after_mm += 1,
                (Del, Del) => del_after_del += 1,
                (In, Del) => del_after_in += 1,
                (Mism, In) | (Match, In) => in_after_mm += 1,
                (In, In) => in_after_in += 1,
                (Del, In) => in_after_del += 1,
                _ => {}
            }
        }
    }
    matchmis -= num_seq;
    let div = |x, y| x as f64 / y as f64;
    //eprintln!("{}\t{}\t{}", matchmis, num_in,num_del);
    Summary::new(
        div(mm_after_mm, matchmis),
        div(in_after_mm, matchmis),
        div(del_after_mm, matchmis),
        div(in_after_in, num_in),
        div(del_after_del, num_del),
        div(del_after_in, num_in),
        div(in_after_del, num_del),
    )
}
