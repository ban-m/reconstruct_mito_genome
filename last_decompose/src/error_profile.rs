use poa_hmm::Config;
use std::collections::HashMap;
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Op {
    Match,
    Mism,
    Del,
    In,
}

pub fn summarize_tab(
    tab: &[last_tiling::LastTAB],
    read: &[bio_utils::fasta::Record],
    contig: &[bio_utils::fasta::Record],
) -> Config {
    let bf = base_freq(&read);
    let reads: HashMap<_, _> = read.iter().map(|e| (e.id().to_string(), e)).collect();
    let contigs: HashMap<_, _> = contig
        .iter()
        .map(|e| (e.id().to_string(), e))
        .collect();
    let ops: Vec<Vec<_>> = tab
        .iter()
        .map(|tab| {
            let rn = tab.seq2_name();
            let cn = tab.seq1_name();
            recover(&reads[rn], &contigs[cn], tab)
        })
        .filter(|ops| ops.len() > 2000)
        .collect();
    summarize_operations(ops, bf)
}

fn base_freq(rs: &[bio_utils::fasta::Record]) -> [f64; 4] {
    let tot = rs.iter().map(|e| e.seq().len()).sum::<usize>() as f64;
    let mut base_count = [0.; 4];
    for base in rs.iter().flat_map(|e| e.seq().iter()) {
        match base {
            b'A' => base_count[0] += 1.,
            b'C' => base_count[1] += 1.,
            b'G' => base_count[2] += 1.,
            b'T' => base_count[3] += 1.,
            _ => {}
        }
    }
    base_count.iter_mut().for_each(|e| *e /= tot);
    base_count
}

fn recover(
    query: &bio_utils::fasta::Record,
    refr: &bio_utils::fasta::Record,
    tab: &last_tiling::LastTAB,
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

fn summarize_operations(opss: Vec<Vec<Op>>, base_freq: [f64; 4]) -> Config {
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
    let _p_in_to_del = div(del_after_in, num_in);
    let p_del_to_in = div(in_after_del, num_del);
    Config {
        mismatch: p_mismatch,
        base_freq,
        p_match,
        p_ins: p_start_in,
        p_del: p_start_del,
        p_extend_ins: p_ext_in,
        p_extend_del: p_ext_del,
        p_del_to_ins: p_del_to_in,
    }
}
