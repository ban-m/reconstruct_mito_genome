#![feature(test)]
extern crate log;
extern crate packed_simd;
extern crate rand;
extern crate rand_xoshiro;
#[cfg(test)]
extern crate rayon;
#[cfg(test)]
extern crate test;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
mod config;
pub use config::*;
pub mod base;
mod base_table;
use base::Base;
mod add_sequence;
mod construct;
mod formatters;
pub mod forward;
pub mod gen_sample;
mod remove_nodes;
const SMALL: f64 = 0.000_000_001;
const LAMBDA: f64 = 0.0;
const FRAC: usize = 10;
const THR: f64 = 1.;
pub mod generate;
#[cfg(test)]
mod tests;
// Edit operation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    Match(usize),
    Deletion(usize),
    Insertion(usize),
    Stop,
}
impl std::fmt::Display for EditOp {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self {
            EditOp::Match(_) => write!(f, "M"),
            EditOp::Deletion(_) => write!(f, "D"),
            EditOp::Insertion(_) => write!(f, "I"),
            EditOp::Stop => write!(f, "S"),
        }
    }
}

pub type POA = PartialOrderAlignment;
#[derive(Clone, Default)]
pub struct PartialOrderAlignment {
    nodes: Vec<Base>,
    weight: f64,
}

type TraceBack = Vec<EditOp>;

impl PartialOrderAlignment {
    pub fn generate(seqs: &[&[u8]], ws: &[f64], config: &Config) -> POA {
        if seqs.is_empty() {
            panic!("Empty string.")
        }
        let seed = (100. * ws.iter().sum::<f64>().floor()) as u64 + seqs.len() as u64;
        let choises: Vec<_> = (0..seqs.len()).collect();
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        let picked = *choises.choose_weighted(&mut rng, |&i| ws[i]).unwrap();
        let (seed, seed_weight) = (&seqs[picked], ws[picked]);
        let max_len = seqs
            .iter()
            .zip(ws.iter())
            .map(|(xs, &w)| if w > 0.001 { xs.len() } else { 0 })
            .max()
            .unwrap_or(0);
        // seqs.iter()
        //     .zip(ws.iter())
        //     .filter(|&(_, &w)| w > 0.001)
        //     .fold(POA::default(), |x, (y, &w)| x.add_with(y, w, config))
        //  .remove_node()
        //     .finalize()
        //     .clean_up()
        seqs.iter()
            .zip(ws.iter())
            .enumerate()
            .filter(|&(idx, (_, &w))| w > 0.001 && idx != picked)
            .map(|(_, x)| x)
            .fold(POA::new(seed, seed_weight), |x, (y, &w)| {
                // .fold(POA::default(), |x, (y, &w)| {
                if x.nodes.len() > 3 * max_len / 2 {
                    x.add_with(y, w, config).remove_node()
                } else {
                    x.add_with(y, w, config)
                }
            })
            .remove_node()
            .clean_up()
            .finalize()
    }
    pub fn update_by(mut self, seqs: &[&[u8]], ws: &[f64], config: &Config) -> POA {
        let max_len = seqs
            .iter()
            .zip(ws.iter())
            .map(|(xs, &w)| if w > 0.001 { xs.len() } else { 0 })
            .max()
            .unwrap_or(0);
        self.weight = 0.;
        self.nodes.iter_mut().for_each(|n| n.weight = 0.);
        seqs.iter()
            .zip(ws.iter())
            .fold(self, |x, (y, &w)| {
                if x.nodes.len() > 3 * max_len / 2 {
                    x.add_with(y, w, config).remove_node()
                } else {
                    x.add_with(y, w, config)
                }
            })
            .remove_node()
            .clean_up()
            .finalize()
    }
    pub fn generate_uniform(seqs: &[&[u8]]) -> POA {
        let ws = vec![1.; seqs.len()];
        POA::generate(seqs, &ws, &DEFAULT_CONFIG)
    }
    pub fn generate_vec(seqs: &[Vec<u8>]) -> POA {
        let ws = vec![1.; seqs.len()];
        let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
        POA::generate(&seqs, &ws, &DEFAULT_CONFIG)
    }
    pub fn new(seq: &[u8], w: f64) -> Self {
        let weight = w;
        let mut nodes: Vec<_> = seq
            .windows(2)
            .enumerate()
            .map(|(idx, pair)| {
                let mut n = Base::new(pair[0]);
                n.add(pair[1], w, idx + 1);
                if idx == 0 {
                    n.is_head = true;
                    n.head_weight += w;
                }
                n
            })
            .collect();
        let mut last = Base::new(*seq.last().unwrap());
        last.is_tail = true;
        nodes.push(last);
        nodes.iter_mut().for_each(|n| n.add_weight(w));
        assert!(nodes.iter().all(|n| n.weight() > 0.));
        Self { weight, nodes }
    }
    //  ----> Graph position ----->
    // 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  |
    //  |
    //Qeury position
    //  |
    //  v
    // Semi-Global alignment containing query sequence.
    pub fn align<F>(&mut self, seq: &[u8], ins: f64, del: f64, score: F) -> (f64, TraceBack)
    where
        F: Fn(u8, u8) -> f64,
    {
        let mut traceback = Vec::with_capacity(seq.len() + 1);
        let edges = self.reverse_edges();
        let mut prev = vec![0.; self.nodes.len() + 1];
        let tb = vec![EditOp::Stop; self.nodes.len() + 1];
        traceback.push(tb);
        let mut updated = vec![];
        let min = -100_000.;
        for (i, &b) in seq.iter().enumerate() {
            updated.push(ins * (i + 1) as f64);
            let mut tb: Vec<_> = Vec::with_capacity(self.nodes.len() + 1);
            tb.push(EditOp::Insertion(0));
            for (j, (n, edges)) in self.nodes.iter().zip(edges.iter()).enumerate() {
                let ms = score(n.base(), b);
                let (m_argmax, m_max, d_argmax, d_max) = if edges.is_empty() {
                    (0, prev[0] + ms, 0, updated[0] + del)
                } else {
                    edges
                        .iter()
                        .map(|&from| from + 1)
                        .map(|idx| (idx, prev[idx] + ms, updated[idx] + del))
                        .fold(
                            (0, min, 0, min),
                            |(m_ax, m_x, d_ax, d_x), (i, m, d)| match (m_x < m, d_x < d) {
                                (true, true) => (i, m, i, d),
                                (true, false) => (i, m, d_ax, d_x),
                                (false, true) => (m_ax, m_x, i, d),
                                (false, false) => (m_ax, m_x, d_ax, d_x),
                            },
                        )
                };
                let max = d_max.max(m_max).max(prev[j + 1] + ins);
                let argmax = match max {
                    x if (x - m_max).abs() < SMALL => EditOp::Match(m_argmax),
                    x if (x - d_max).abs() < SMALL => EditOp::Deletion(d_argmax),
                    _ => EditOp::Insertion(j + 1),
                };
                // Select one of the optimal operations.
                updated.push(max);
                tb.push(argmax)
            }
            traceback.push(tb);
            std::mem::swap(&mut prev, &mut updated);
            updated.clear();
        }
        // eprintln!("OK");
        // Traceback
        // g_pos = position on the graph, q_pos = position on the query
        let mut q_pos = seq.len();
        let mut operations = vec![];
        let (mut g_pos, poa_score) = prev
            .into_iter()
            .enumerate()
            .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
            .unwrap_or_else(|| panic!("{}", line!()));
        while traceback[q_pos][g_pos] != EditOp::Stop {
            match traceback[q_pos][g_pos] {
                EditOp::Match(from) => {
                    operations.push(EditOp::Match(g_pos - 1));
                    g_pos = from;
                    q_pos -= 1;
                }
                EditOp::Deletion(from) => {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = from;
                }
                EditOp::Insertion(_) => {
                    operations.push(EditOp::Insertion(g_pos - 1));
                    q_pos -= 1;
                }
                EditOp::Stop => break,
            }
        }
        operations.reverse();
        (poa_score, operations)
    }
    pub fn add_with(self, seq: &[u8], w: f64, c: &Config) -> Self {
        let ins = c.p_ins.ln();
        let del = c.p_del.ln();
        let mat = -10. * c.p_match.ln();
        let mism = c.mismatch.ln();
        self.add_w_param(seq, w, ins, del, |x, y| if x == y { mat } else { mism })
    }
    pub fn add(self, seq: &[u8], w: f64) -> Self {
        self.add_w_param(seq, w, -1., -1., |x, y| if x == y { 1. } else { -1. })
    }
    pub fn add_w_param<F>(mut self, seq: &[u8], w: f64, ins: f64, del: f64, score: F) -> Self
    where
        F: Fn(u8, u8) -> f64,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, w);
        }
        // Alignment
        let (_, traceback) = self.align(seq, ins, del, score);
        let mut q_pos = 0;
        let mut previous: Option<usize> = None;
        for operation in traceback {
            match operation {
                EditOp::Match(to) => {
                    let base = seq[q_pos];
                    let position = if self.nodes[to].base() == base {
                        to
                    } else {
                        let mut new_node = Base::new(seq[q_pos]);
                        new_node.ties.push(to);
                        let position = self.nodes.len();
                        self.nodes.push(new_node);
                        self.nodes[to].ties.push(position);
                        position
                    };
                    self.nodes[position].add_weight(w);
                    if q_pos == seq.len() - 1 {
                        self.nodes[position].is_tail = true;
                        self.nodes[position].tail_weight += w;
                    } else if q_pos == 0 {
                        self.nodes[position].is_head = true;
                        self.nodes[position].head_weight += w;
                    }
                    if let Some(p) = previous {
                        self.nodes[p].add(base, w, position);
                    };
                    previous = Some(position);
                    q_pos += 1;
                }
                EditOp::Insertion(_) => {
                    let base = seq[q_pos];
                    let mut new_node = Base::new(base);
                    new_node.add_weight(w);
                    if q_pos == seq.len() - 1 {
                        new_node.is_tail = true;
                        new_node.tail_weight += w;
                    } else if q_pos == 0 {
                        new_node.is_head = true;
                        new_node.head_weight += w;
                    }
                    self.nodes.push(new_node);
                    if let Some(p) = previous {
                        let position = self.nodes.len() - 1;
                        self.nodes[p].add(base, w, position);
                    }
                    previous = Some(self.nodes.len() - 1);
                    q_pos += 1;
                }
                EditOp::Deletion(_) | EditOp::Stop => {}
            }
        }
        assert_eq!(q_pos, seq.len());
        self.weight += w;
        assert!(self.nodes.iter().all(|node| node.weight() > 0.));
        self.topological_sort()
    }
    // Return the minimum distance from root node to each node.
    // Done by BFS.
    // pub fn min_dist(&self) -> Vec<(usize, usize)> {
    //     let mut dist = vec![(0, 0); self.nodes.len()];
    //     let mut queue: Vec<_> = {
    //         let mut is_head = vec![true; self.nodes.len()];
    //         for &to in self.nodes.iter().flat_map(|n| n.edges.iter()) {
    //             is_head[to] = false;
    //         }
    //         is_head
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(idx, &e)| if e { Some(idx) } else { None })
    //             .collect()
    //     };
    //     let mut flag: Vec<_> = vec![false; self.nodes.len()];
    //     let mut depth = 0;
    //     let mut update = vec![];
    //     while !queue.is_empty() {
    //         while let Some(from) = queue.pop() {
    //             if flag[from] {
    //                 continue;
    //             }
    //             flag[from] = true;
    //             dist[from].0 = depth;
    //             for &to in self.nodes[from].edges.iter() {
    //                 if !flag[to] {
    //                     dist[to].1 = from;
    //                 }
    //                 update.push(to);
    //             }
    //         }
    //         std::mem::swap(&mut queue, &mut update);
    //         depth += 1;
    //     }
    //     dist
    // }
    pub fn edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[from].push(to);
            }
            for &tied in n.ties.iter() {
                for &to in self.nodes[tied].edges.iter() {
                    edges[from].push(to);
                }
            }
        }
        edges
    }
    pub fn reverse_edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[to].push(from);
            }
            for &tied in n.ties.iter() {
                for &to in self.nodes[tied].edges.iter() {
                    edges[to].push(from);
                }
            }
        }
        edges
    }
    pub fn clean_up(self) -> Self {
        self
    }
    pub fn finalize(mut self) -> Self {
        self.nodes.iter_mut().for_each(|e| e.finalize());
        self
    }
}
