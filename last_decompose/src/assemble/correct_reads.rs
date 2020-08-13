use rayon::prelude::*;
use std::collections::HashMap;
fn score((q_u, q_c): (u64, u64), (r_u, r_c): (u64, u64), mat: i32, mism: i32) -> i32 {
    if q_u == r_u && q_c == r_c {
        mat
    } else if q_u == r_u {
        mism
    } else {
        -100
    }
}

// Align the query to the reference and
// return the edit operations. Note that
// A Cigar::Match(x,y) mean the query sequence at that point is (x,y)
// And Cigar::Ins is a insertion to the reference.
// Also, the alignment is "semi-global" one. See the initialization step.
// TODO: faster!
fn alignment(
    qry: &[(u64, u64)],
    rfr: &[(u64, u64)],
    (mat, mism, gap): (i32, i32, i32),
) -> Option<(i32, Vec<Cigar>)> {
    let mut dp = vec![vec![0; rfr.len() + 1]; qry.len() + 1];
    for (i, &q) in qry.iter().enumerate() {
        for (j, &r) in rfr.iter().enumerate() {
            let mat = dp[i][j] + score(q, r, mat, mism);
            let ins = dp[i][j + 1] + gap;
            let del = dp[i + 1][j] + gap;
            dp[i + 1][j + 1] = mat.max(ins).max(del);
        }
    }
    // Determine the starting point.
    let (row_pos, row_max) = dp.last()?.iter().enumerate().max_by_key(|x| x.1)?;
    let (column_pos, column_max) = dp
        .iter()
        .filter_map(|x| x.last())
        .enumerate()
        .max_by_key(|x| x.1)?;
    let score = *column_max.max(row_max);
    if score <= mat {
        return None;
    }
    let (mut q_pos, mut r_pos) = if row_max < column_max {
        (column_pos, rfr.len())
    } else {
        (qry.len(), row_pos)
    };
    assert_eq!(dp[q_pos][r_pos], *column_max.max(row_max));
    let mut cigar = vec![];
    for q in (q_pos..qry.len()).rev() {
        let (unit, cluster) = qry[q];
        cigar.push(Cigar::Ins(unit, cluster));
    }
    for _ in (r_pos..rfr.len()).rev() {
        cigar.push(Cigar::Del);
    }
    // Traceback.
    while q_pos > 0 && r_pos > 0 {
        let current = dp[q_pos][r_pos];
        let op = if current == dp[q_pos - 1][r_pos] + gap {
            let (unit, cluster) = qry[q_pos - 1];
            q_pos -= 1;
            Cigar::Ins(unit, cluster)
        } else if current == dp[q_pos][r_pos - 1] + gap {
            r_pos -= 1;
            Cigar::Del
        } else {
            let (unit, cluster) = qry[q_pos - 1];
            r_pos -= 1;
            q_pos -= 1;
            Cigar::Match(unit, cluster)
        };
        cigar.push(op);
    }
    while q_pos > 0 {
        let (unit, cluster) = qry[q_pos - 1];
        cigar.push(Cigar::Ins(unit, cluster));
        q_pos -= 1;
    }
    while r_pos > 0 {
        cigar.push(Cigar::Del);
        r_pos -= 1;
    }
    cigar.reverse();
    Some((score, cigar))
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Cigar {
    Match(u64, u64),
    Ins(u64, u64),
    Del,
}

#[derive(Clone, Debug)]
struct Pileup {
    column: Vec<Column>,
}

impl Pileup {
    fn new(ns: &[(u64, u64)]) -> Self {
        let column: Vec<_> = ns.iter().map(|&n| Column::new(n)).collect();
        Self { column }
    }
    fn add(mut self, y: Vec<Cigar>) -> Self {
        let mut r_pos = 0;
        for op in y {
            match op {
                Cigar::Match(unit, cl) => {
                    self.column[r_pos].m.push((unit, cl));
                    r_pos += 1;
                }
                Cigar::Ins(_, _) => {}
                Cigar::Del => {
                    r_pos += 1;
                }
            }
        }
        self
    }
}

#[derive(Clone, Debug)]
struct Column {
    m: Vec<(u64, u64)>,
}

impl Column {
    fn new(n: (u64, u64)) -> Self {
        Self { m: vec![n] }
    }
    fn generate(&self) -> (u64, u64) {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for &x in self.m.iter() {
            *counts.entry(x).or_default() += 1;
        }
        counts.into_iter().max_by_key(|x| x.1).unwrap().0
    }
}

pub fn correct_reads(reads: &mut [super::ChunkedRead]) {
    let reads_summary: Vec<Vec<(u64, u64)>> = reads
        .iter()
        .map(|read| {
            read.nodes
                .iter()
                .map(|n| (n.window_position as u64, n.cluster as u64))
                .collect()
        })
        .collect();
    let rev_for_reads: Vec<_> = {
        let rev = reads_summary
            .iter()
            .map(|read| read.iter().copied().rev().collect::<Vec<_>>());
        reads_summary.iter().cloned().zip(rev).collect()
    };
    let corrected_reads: Vec<_> = reads_summary
        .par_iter()
        .map(|read| correct_read(read, &rev_for_reads))
        .collect();
    assert_eq!(reads.len(), corrected_reads.len());
    for (read, corrected) in reads.iter_mut().zip(corrected_reads) {
        assert_eq!(read.nodes.len(), corrected.len());
        for (node, (unit, cluster)) in read.nodes.iter_mut().zip(corrected) {
            node.window_position = unit as usize;
            node.cluster = cluster as u8;
        }
    }
}

fn correct_read(
    read: &[(u64, u64)],
    reads: &[(Vec<(u64, u64)>, Vec<(u64, u64)>)],
) -> Vec<(u64, u64)> {
    let param = (1, -1, -3);
    let pileup = reads
        .iter()
        .filter_map(|(forward, rev)| match alignment(forward, read, param) {
            Some(res) => Some(res),
            None => alignment(rev, read, param),
        })
        .fold(Pileup::new(read), |x, (_, y)| x.add(y));
    pileup
        .column
        .into_iter()
        .map(|column| column.generate())
        .collect()
}
