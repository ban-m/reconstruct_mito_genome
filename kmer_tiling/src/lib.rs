extern crate rand;
#[macro_use]
extern crate log;
extern crate env_logger;
use rand::{
    distributions::{Bernoulli, Distribution},
    thread_rng,
};
const BASE: &[u8] = b"ACGT";

type Kmer = Vec<u8>;
/// Enumerate all kmers in distance `d` from the given kmer.

#[derive(PartialEq, Eq, Clone)]
enum Op {
    Del,
    Ins,
    Subst,
    Match,
}

impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Op::Del => write!(f, "X"),
            Op::Ins => write!(f, "V"),
            Op::Subst => write!(f, "S"),
            Op::Match => write!(f, "O"),
        }
    }
}
impl std::fmt::Debug for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

#[allow(dead_code)]
fn to_string(op: &[Op]) -> String {
    op.iter()
        .map(|op| match op {
            &Op::Ins => 'V',
            &Op::Del => 'X',
            &Op::Subst => 'S',
            &Op::Match => 'O',
        })
        .collect()
}

pub fn enumerate_dist_k(kmer: &[u8], d: u16) -> Vec<Kmer> {
    let max_indel = d / 2;
    info!(
        "Enumerate all kmer within dist {} from {}",
        d,
        String::from_utf8_lossy(kmer)
    );
    info!("max indel:{}", max_indel);
    (0..=max_indel)
        .flat_map(|indel| enumerate_kmer_with_indel(kmer, indel, d - 2 * indel))
        .filter(|cand| edit_dist(cand, kmer) == d)
        .collect()
}

pub fn enumerate_dist_k_squished(kmer: &[u8], d: u16) -> Vec<Kmer> {
    // Indel would be treated as half penalty as mismatch.
    let max_indel = d;
    (0..=max_indel)
        .flat_map(|indel| enumerate_kmer_with_indel(kmer, indel, d - indel))
        .filter(|cand| edit_dist_with(cand, kmer, 2, 1) == 2 * d)
        .collect()
}

pub fn enumerate_dist_k_squished_withprob(kmer: &[u8], d: u16, p: f64) -> Vec<Kmer> {
    let max_indel = d;
    (0..=max_indel)
        .flat_map(|indel| enumerate_kmer_with_indel_withprob(kmer, indel, d - indel, p))
        .filter(|cand| edit_dist_with(cand, kmer, 2, 1) == 2 * d)
        .collect()
}

fn enumerate_kmer_with_indel(kmer: &[u8], indel: u16, subst: u16) -> Vec<Kmer> {
    // the number of insertion and delpetion should be the same.
    // The patterns should be the combination of positions and bases.
    // Notice that, there should not be any insertion followed by deletion or
    // deletion followed by insertion, as it is better to treat as one single substitution.
    let (del, ins, sub, mat) = (indel, indel, subst, kmer.len() as u16 - indel - subst);
    enumerate_all_edit_pattern(del, ins, sub, mat)
        .into_iter()
        .flat_map(|pat| enumerate_modifications(kmer, pat))
        .collect()
}

fn enumerate_kmer_with_indel_withprob(kmer: &[u8], indel: u16, subst: u16, p: f64) -> Vec<Kmer> {
    // the number of insertion and delpetion should be the same.
    // The patterns should be the combination of positions and bases.
    // Notice that, there should not be any insertion followed by deletion or
    // deletion followed by insertion, as it is better to treat as one single substitution.
    let (del, ins, sub, mat) = (indel, indel, subst, kmer.len() as u16 - indel - subst);
    let mut rng = thread_rng();
    let ber = Bernoulli::new(p);
    enumerate_all_edit_pattern(del, ins, sub, mat)
        .into_iter()
        .filter(|_| ber.sample(&mut rng))
        .flat_map(|pat| {
            enumerate_modifications(kmer, pat)
        })
        .collect()
}

fn enumerate_all_edit_pattern(del: u16, ins: u16, sub: u16, mat: u16) -> Vec<Vec<Op>> {
    enumerate_all_edit_pattern_inner(del, ins, sub, mat, vec![])
}

fn enumerate_all_edit_pattern_inner(
    del: u16,
    ins: u16,
    sub: u16,
    mat: u16,
    current: Vec<Op>,
) -> Vec<Vec<Op>> {
    // info!(
    //     "Ops:{},(d,i,s,m) = ({},{},{},{})",
    //     to_string(&current),
    //     del,
    //     ins,
    //     sub,
    //     mat
    // );
    if del == 0 && ins == 0 && sub == 0 && mat == 0 {
        vec![current]
    } else {
        let mut res = vec![];
        // Branching.
        // The match and substitution is "context-free", in that it can be successor of any other oparations.
        if sub != 0 {
            let mut current = current.clone();
            current.push(Op::Subst);
            res.extend(enumerate_all_edit_pattern_inner(
                del,
                ins,
                sub - 1,
                mat,
                current,
            ));
        }
        if mat != 0 {
            let mut current = current.clone();
            current.push(Op::Match);
            res.extend(enumerate_all_edit_pattern_inner(
                del,
                ins,
                sub,
                mat - 1,
                current,
            ));
        }
        // Remark that Del-Ins and Ins-Del are forbidden.
        let (del_ok, ins_ok) = if current.len() == 0 {
            (del != 0, ins != 0)
        } else {
            match current.last().unwrap() {
                &Op::Del => (del != 0, false), // successive del is ok.
                &Op::Ins => (false, ins != 0),
                _ => (del != 0, ins != 0),
            }
        };
        if del_ok {
            let mut current = current.clone();
            current.push(Op::Del);
            res.extend(enumerate_all_edit_pattern_inner(
                del - 1,
                ins,
                sub,
                mat,
                current,
            ));
        }
        if ins_ok {
            let mut current = current.clone();
            current.push(Op::Ins);
            res.extend(enumerate_all_edit_pattern_inner(
                del,
                ins - 1,
                sub,
                mat,
                current,
            ));
        }
        res
    }
}

fn enumerate_modifications(kmer: &[u8], pat: Vec<Op>) -> Vec<Kmer> {
    // Intoroduce all possible modifications specified by pat on the given kmer.
    enumerate_modifications_inner(kmer, &pat, vec![])
}

fn enumerate_modifications_inner(kmer: &[u8], pat: &[Op], current: Vec<u8>) -> Vec<Kmer> {
    if pat.is_empty() {
        vec![current]
    } else {
        let mut res = vec![];
        let (op, pat) = pat.split_first().unwrap();
        if op == &Op::Del {
            res.extend(enumerate_modifications_inner(&kmer[1..], pat, current));
        } else if op == &Op::Ins {
            for base in BASE {
                let mut current = current.clone();
                current.push(*base);
                res.extend(enumerate_modifications_inner(kmer, pat, current));
            }
        } else if op == &Op::Subst {
            let from = kmer[0];
            for base in BASE {
                if base != &from {
                    let mut current = current.clone();
                    current.push(*base);
                    res.extend(enumerate_modifications_inner(&kmer[1..], pat, current));
                }
            }
        } else if op == &Op::Match {
            let mut current = current.clone();
            current.push(kmer[0]);
            res.extend(enumerate_modifications_inner(&kmer[1..], pat, current));
        }
        res
    }
}

pub fn edit_dist(x: &[u8], y: &[u8]) -> u16 {
    let mut d = vec![vec![0; y.len() + 1]; x.len() + 1];
    for i in 0..=x.len() {
        d[i][0] = i as u16;
    }
    for j in 0..=y.len() {
        d[0][j] = j as u16;
    }
    for i in 1..=x.len() {
        for j in 1..=y.len() {
            let is_match = if x[i - 1] == y[j - 1] { 0 } else { 1 };
            d[i][j] = (d[i - 1][j - 1] + is_match)
                .min(d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1);
        }
    }
    d[x.len()][y.len()]
}

pub fn edit_dist_with(x: &[u8], y: &[u8], mis: u16, indel: u16) -> u16 {
    let mut d = vec![vec![0; y.len() + 1]; x.len() + 1];
    for i in 0..=x.len() {
        d[i][0] = i as u16;
    }
    for j in 0..=y.len() {
        d[0][j] = j as u16;
    }
    for i in 1..=x.len() {
        for j in 1..=y.len() {
            let is_match = if x[i - 1] == y[j - 1] { 0 } else { mis };
            d[i][j] = (d[i - 1][j - 1] + is_match)
                .min(d[i - 1][j] + indel)
                .min(d[i][j - 1] + indel);
        }
    }
    d[x.len()][y.len()]
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn enumerate_dist_k_test() {
        let kmer = b"AAAAA".to_vec();
        let kmers = enumerate_dist_k(&kmer, 0);
        assert_eq!(vec![kmer], kmers);
        let kmer = b"CTGC".to_vec();
        let kmers = enumerate_dist_k(&kmer, 0);
        assert_eq!(vec![kmer], kmers);
        let kmer = b"AAAA".to_vec();
        let mut kmers = enumerate_dist_k(&kmer, 1);
        kmers.sort();
        let answer: Vec<_> = vec![
            b"AAAC", b"AAAG", b"AAAT", b"AACA", b"AAGA", b"AATA", b"ACAA", b"AGAA", b"ATAA",
            b"CAAA", b"GAAA", b"TAAA",
        ]
        .into_iter()
        .map(|e| e.to_vec())
        .collect();;
        assert_eq!(kmers, answer);
        let kmers = enumerate_dist_k(&kmer, 2);
        for res in kmers {
            debug_assert!(
                edit_dist(&kmer, &res) <= 2,
                "{}-{}",
                String::from_utf8_lossy(&kmer),
                String::from_utf8_lossy(&res)
            );
        }
        let kmer = b"ACATCAGT";
        for res in enumerate_dist_k(kmer, 2) {
            assert!(edit_dist(kmer, &res) <= 2);
        }
        for res in enumerate_dist_k(kmer, 5) {
            assert!(edit_dist(kmer, &res) <= 5);
        }
    }
    #[test]
    fn enumerate_pat_test() {
        use Op::*;
        assert_eq!(
            vec![vec![Op::Match; 5]],
            enumerate_all_edit_pattern(0, 0, 0, 5)
        );
        assert_eq!(
            vec![vec![Op::Del; 4]],
            enumerate_all_edit_pattern(4, 0, 0, 0)
        );
        assert_eq!(
            vec![vec![Op::Ins; 4]],
            enumerate_all_edit_pattern(0, 4, 0, 0)
        );
        assert_eq!(
            vec![vec![Op::Subst; 4]],
            enumerate_all_edit_pattern(0, 0, 4, 0)
        );
        let res = enumerate_all_edit_pattern(1, 0, 0, 4);
        for test in vec![
            vec![Match, Match, Match, Match, Del],
            vec![Match, Match, Match, Del, Match],
            vec![Match, Match, Del, Match, Match],
            vec![Match, Del, Match, Match, Match],
            vec![Del, Match, Match, Match, Match],
        ] {
            assert!(res.contains(&test));
        }
        let res = enumerate_all_edit_pattern(0, 0, 1, 3);
        for test in vec![
            vec![Match, Match, Match, Subst],
            vec![Match, Match, Subst, Match],
            vec![Match, Subst, Match, Match],
            vec![Subst, Match, Match, Match],
        ] {
            assert!(res.contains(&test));
        }
        let res = enumerate_all_edit_pattern(1, 1, 0, 3);
        for test in vec![
            vec![Match, Match, Ins, Match, Del],
            vec![Match, Ins, Match, Match, Del],
            vec![Ins, Match, Match, Match, Del],
            vec![Match, Ins, Match, Del, Match],
            vec![Ins, Match, Match, Del, Match],
            vec![Ins, Match, Del, Match, Match],
            vec![Match, Match, Del, Match, Ins],
            vec![Match, Del, Match, Match, Ins],
            vec![Match, Del, Match, Ins, Match],
            vec![Del, Match, Match, Match, Ins],
            vec![Del, Match, Match, Ins, Match],
            vec![Del, Match, Ins, Match, Match],
        ] {
            assert!(res.contains(&test));
        }
        for neg_test in vec![
            vec![Match, Match, Match, Ins, Del],
            vec![Match, Match, Match, Del, Ins],
            vec![Match, Match, Ins, Del, Match],
            vec![Match, Ins, Del, Match, Match],
            vec![Match, Del, Ins, Match, Match],
        ] {
            assert!(!res.contains(&neg_test));
        }
    }
}
