extern crate bio_utils;
extern crate last_tiling;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let tab: Vec<_> = last_tiling::parse_tab_file(&args[1])?;
    let read: Vec<_> = bio_utils::fasta::parse_into_vec(&args[2])?;
    let contig: Vec<_> = bio_utils::fasta::parse_into_vec(&args[3])?;
    eprintln!("{}", tab.len());
    {
        let summary = summarize_tab(tab, read, contig);
        println!("{}", summary);
    }
    Ok(())
}

// TODO: write function to convert Sumamry to Config.
#[derive(Debug, Clone, Copy)]
struct Summary {
    // Note: this is a mismatch /  match + mismatch.
    p_mismatch: f64,
    p_match: f64,
    p_start_in: f64,
    p_start_del: f64,
    p_ext_in: f64,
    p_ext_del: f64,
    p_in_to_del: f64,
    p_del_to_in: f64,
    base_freq: [f64; 4],
}

impl Summary {}

impl std::fmt::Display for Summary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Pmism\tPm\tPSi\tPSd\tPEi\tPEd\tPid\tPdi")?;
        writeln!(
            f,
            "{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
            100. * self.p_mismatch,
            100. * self.p_match,
            100. * self.p_start_in,
            100. * self.p_start_del,
            100. * self.p_ext_in,
            100. * self.p_ext_del,
            100. * self.p_in_to_del,
            100. * self.p_del_to_in
        )?;
        writeln!(f, "A\tC\tG\tT")?;
        write!(
            f,
            "{:.2}\t{:.2}\t{:.2}\t{:.2}",
            100. * self.base_freq[0],
            100. * self.base_freq[1],
            100. * self.base_freq[2],
            100. * self.base_freq[3]
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

fn base_freq(rs: &[bio_utils::fasta::Record]) -> [f64; 4] {
    let tot = rs.iter().map(|e| e.seq().len()).sum::<usize>() as f64;
    let bf = rs
        .iter()
        .map(|e| {
            e.seq().iter().fold([0; 4], |mut x, b| {
                match b {
                    b'A' => x[0] += 1,
                    b'C' => x[1] += 1,
                    b'G' => x[2] += 1,
                    b'T' => x[3] += 1,
                    _ => {}
                };
                x
            })
        })
        .fold([0; 4], |mut x, e| {
            x[0] += e[0];
            x[1] += e[1];
            x[2] += e[2];
            x[3] += e[3];
            x
        });
    [
        bf[0] as f64 / tot,
        bf[1] as f64 / tot,
        bf[2] as f64 / tot,
        bf[3] as f64 / tot,
    ]
}

fn summarize_tab(
    tab: Vec<last_tiling::LastTAB>,
    read: Vec<bio_utils::fasta::Record>,
    contig: Vec<bio_utils::fasta::Record>,
) -> Summary {
    let bf = base_freq(&read);
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
        .filter(|ops| ops.len() > 2000)
        .collect();
    eprintln!("{}", ops.len());
    summarize_operations(ops, bf)
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

fn summarize_operations(opss: Vec<Vec<Op>>, base_freq: [f64; 4]) -> Summary {
    // match + mismatch.
    let mut matchmis = 0;
    let mut num_mis = 0;
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
        num_mis += ops
            .iter()
            .filter(|&e| match e {
                Mism => true,
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
    let div = |x, y| x as f64 / y as f64;
    let p_mismatch = div(num_mis, matchmis);
    matchmis -= num_seq;
    //eprintln!("{}\t{}\t{}", matchmis, num_in,num_del);
    let p_match = div(mm_after_mm, matchmis);
    let p_start_in = div(in_after_mm, matchmis);
    let p_start_del = div(del_after_mm, matchmis);
    let p_ext_in = div(in_after_in, num_in);
    let p_ext_del = div(del_after_del, num_del);
    let p_in_to_del = div(del_after_in, num_in);
    let p_del_to_in = div(in_after_del, num_del);
    Summary {
        p_mismatch,
        p_match,
        p_start_in,
        p_start_del,
        p_ext_in,
        p_ext_del,
        p_in_to_del,
        p_del_to_in,
        base_freq,
    }
}
