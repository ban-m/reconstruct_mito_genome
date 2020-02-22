#![feature(test)]
extern crate log;
extern crate packed_simd;
extern crate rand;
extern crate rand_xoshiro;
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
pub mod generate;
#[cfg(test)]
mod tests;
// Edit operation
#[derive(Debug, Clone, Copy)]
pub enum EditOp {
    Match(usize),
    Deletion(usize),
    Insertion(usize),
    Stop,
}

pub type POA = PartialOrderAlignment;
#[derive(Clone, Default)]
pub struct PartialOrderAlignment {
    nodes: Vec<Base>,
    weight: f64,
}

type DP = Vec<f64>;
type TraceBack = Vec<Vec<EditOp>>;

impl PartialOrderAlignment {
    pub fn generate(seqs: &[&[u8]], ws: &[f64], config: &Config) -> POA {
        if seqs.is_empty() {
            panic!("Empty string.")
        }
        let seed = ws.iter().sum::<f64>().floor() as u64 + seqs.len() as u64;
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
        seqs.iter()
            .zip(ws.iter())
            .fold(POA::default(), |x, (y, &w)| x.add_with(y, w, config))
        // seqs.iter()
        //     .zip(ws.iter())
        //     .enumerate()
        //     .filter(|&(idx, (_, &w))| w > 0.001 && idx != picked)
        //     .map(|(_, x)| x)
        //     .fold((POA::new(seed, seed_weight), 0), |(x, t), (y, &w)| {
        //         if (t + 1) % 20 != 0 || t + 20 > seqs.len() {
        //             (x.add_with(y, w, config), t + 1)
        //         } else {
        //             let x = x
        //                 .add_with(y, w, config)
        //                 .remove_node(max_len)
        //                 .select_head_tail();
        //             (x, t + 1)
        //         }
        //     })
        //     .0
        //     // .remove_node(max_len)
        //     // .select_head_tail()
        //     .clean_up()
        //     .finalize()
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
    //  |
    //  |
    //Qeury position
    //  |
    //  v
    pub fn align<F>(&mut self, seq: &[u8], ins: f64, del: f64, score: F) -> (DP, TraceBack)
    where
        F: Fn(u8, u8) -> f64,
    {
        let mut traceback = vec![];
        let edges = self.reverse_edges();
        assert!(self.nodes.iter().any(|e| e.is_head));
        let (mut prev, mut tb): (Vec<_>, Vec<_>) = self
            .min_dist()
            .into_iter()
            .map(|(dist, parent)| (del * (dist + 1) as f64, EditOp::Deletion(parent)))
            .unzip();
        prev.insert(0, 0.);
        tb.insert(0, EditOp::Stop);
        traceback.push(tb);
        let mut updated = vec![];
        let min = -100_000.;
        for (i, &b) in seq.iter().enumerate() {
            updated.push(ins * (i + 1) as f64);
            let mut tb = Vec::with_capacity(seq.len() + 1);
            tb.push(EditOp::Insertion(0));
            tb.extend(self.nodes.iter().enumerate().map(|(j, n)| {
                let ms = score(n.base(), b);
                let (m_argmax, m_max, d_argmax, d_max) = if edges[j].is_empty() {
                    (0, prev[0] + ms, 0, updated[0] + del)
                } else if edges[j].len() == 1 {
                    let u = edges[j][0];
                    (u + 1, prev[u + 1] + ms, u + 1, updated[u + 1] + del)
                } else {
                    edges[j]
                        .iter()
                        .map(|&from| (from + 1, prev[from + 1] + ms, updated[from + 1] + del))
                        .fold(
                            (0, min, 0, min),
                            |(mut m_ax, mut m_x, mut d_ax, mut d_x), (i, m, d)| {
                                if m_x < m {
                                    m_ax = i;
                                    m_x = m;
                                }
                                if d_x < d {
                                    d_ax = i;
                                    d_x = d;
                                }
                                (m_ax, m_x, d_ax, d_x)
                            },
                        )
                };
                // Insertion into the reference(graph)
                let (i_argmax, i_max) = (j + 1, prev[j + 1] + ins);
                // Select one of the optimal operations.
                let (argmax, max) = if d_max <= m_max && i_max <= m_max {
                    (EditOp::Match(m_argmax), m_max)
                } else if m_max <= d_max && i_max <= d_max {
                    (EditOp::Deletion(d_argmax), d_max)
                } else {
                    (EditOp::Insertion(i_argmax), i_max)
                };
                // Fill.
                updated.push(max);
                argmax
            }));
            traceback.push(tb);
            std::mem::swap(&mut prev, &mut updated);
            updated.clear();
        }
        (prev, traceback)
    }
    pub fn add_with(self, seq: &[u8], w: f64, c: &Config) -> Self {
        let ins = (c.p_ins + c.p_del).ln();
        let del = ins;
        let mism = (c.mismatch + c.mismatch).ln();
        let mat = (c.p_match - c.p_ins - c.p_del - c.mismatch).ln();
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
        let (dp, traceback) = self.align(seq, ins, del, score);
        // Traceback
        // g_pos = position on the graph, q_pos = position on the query
        assert_eq!(self.nodes.len() + 1, dp.len());
        let mut q_pos = seq.len();
        let (mut g_pos, _) = self
            .nodes
            .iter()
            .zip(dp.into_iter().skip(1))
            .enumerate()
            .flat_map(|(idx, (node, score))| {
                if node.is_tail {
                    Some((idx + 1, score))
                } else {
                    None
                }
            })
            .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
            .unwrap();
        let mut previous = None;
        while q_pos > 0 || g_pos > 0 {
            eprintln!("{:?}", traceback[q_pos][g_pos]);
            match traceback[q_pos][g_pos] {
                EditOp::Match(from) => {
                    let current_pos = if self.nodes[g_pos - 1].base() == seq[q_pos - 1] {
                        g_pos - 1
                    } else {
                        self.nodes.push(Base::new(seq[q_pos - 1]));
                        self.nodes.len() - 1
                    };
                    self.nodes[current_pos].add_weight(w);
                    if q_pos == seq.len() {
                        self.nodes[current_pos].is_tail = true;
                    } else if q_pos == 1 {
                        self.nodes[current_pos].is_head = true;
                    }
                    if let Some(p) = previous {
                        let base = self.nodes[p as usize].base();
                        self.nodes[current_pos].add(base, w, p);
                    };
                    previous = Some(current_pos);
                    g_pos = from;
                    q_pos -= 1;
                }
                EditOp::Deletion(from) => {
                    // We do not increment the weight of g_pos,
                    // because the weight is query's, not graph's.
                    g_pos = from;
                }
                EditOp::Insertion(_) => {
                    let mut new_node = Base::new(seq[q_pos - 1]);
                    new_node.add_weight(w);
                    if q_pos == seq.len() {
                        new_node.is_tail = true;
                    } else if q_pos == 1 {
                        new_node.is_head = true;
                    }
                    self.nodes.push(new_node);
                    if let Some(p) = previous {
                        let base = self.nodes[p].base();
                        self.nodes.last_mut().unwrap().add(base, w, p);
                    }
                    previous = Some(self.nodes.len() - 1);
                    q_pos -= 1;
                }
                EditOp::Stop => break,
            }
        }
        self.weight += w;
        assert!(self.nodes.iter().all(|node| node.weight() > 0.));
        self.topological_sort()
    }
    // Return the minimum distance from root node to each node.
    // Done by BFS.
    pub fn min_dist(&self) -> Vec<(usize, usize)> {
        let mut dist = vec![(0, 0); self.nodes.len()];
        let mut queue: Vec<_> = {
            let mut is_head = vec![true; self.nodes.len()];
            for &to in self.nodes.iter().flat_map(|n| n.edges.iter()) {
                is_head[to] = false;
            }
            is_head
                .iter()
                .enumerate()
                .filter_map(|(idx, &e)| if e { Some(idx) } else { None })
                .collect()
        };
        let mut flag: Vec<_> = vec![false; self.nodes.len()];
        let mut depth = 0;
        let mut update = vec![];
        while !queue.is_empty() {
            while let Some(from) = queue.pop() {
                if flag[from] {
                    continue;
                }
                flag[from] = true;
                dist[from].0 = depth;
                for &to in self.nodes[from].edges.iter() {
                    if !flag[to] {
                        dist[to].1 = from;
                    }
                    update.push(to);
                }
            }
            std::mem::swap(&mut queue, &mut update);
            depth += 1;
        }
        dist
    }
    pub fn edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[from].push(to);
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
    fn select_head_tail(mut self) -> Self {
        let mut is_head = vec![true; self.nodes.len()];
        for n in self.nodes.iter() {
            for &to in n.edges.iter() {
                is_head[to] = false;
            }
        }
        assert!(is_head.iter().any(|&e| e), "{:?}", self);
        if self.nodes.iter().all(|e| !e.is_head) {
            self.nodes.iter_mut().zip(is_head).for_each(|(n, is_head)| {
                n.is_head = is_head;
            });
        }
        if self.nodes.iter().all(|e| !e.is_tail) {
            self.nodes.iter_mut().for_each(|n| {
                n.is_tail = n.edges.is_empty();
            });
        }
        assert!(self.nodes.iter().any(|n| n.is_head));
        assert!(self.nodes.iter().any(|n| n.is_tail));
        self
    }
}
