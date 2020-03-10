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
const LAMBDA_INS: f64 = 0.05;
const LAMBDA_MATCH: f64 = 0.1;
const THR: f64 = 0.4;
const MIN: i32 = -100000;
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
    dp: Vec<Vec<i32>>,
    route_weight: Vec<Vec<f32>>,
    last_elements: Vec<i32>,
}

type TraceBack = Vec<EditOp>;

impl PartialOrderAlignment {
    pub fn nodes(&self) -> &[Base] {
        &self.nodes
    }
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    pub fn num_edges(&self) -> usize {
        self.nodes.iter().map(|n| n.edges.len()).sum::<usize>()
    }
    pub fn view(&self, seq: &[u8], traceback: &[EditOp]) -> (String, String) {
        let mut q_pos = 0;
        let (mut q, mut g) = (String::new(), String::new());
        for &op in traceback {
            match op {
                EditOp::Deletion(g_pos) => {
                    q.push('-');
                    g.push(self.nodes()[g_pos].base() as char);
                }
                EditOp::Insertion(_) => {
                    g.push('-');
                    q.push(seq[q_pos] as char);
                    q_pos += 1;
                }
                EditOp::Match(g_pos) => {
                    g.push(self.nodes()[g_pos].base() as char);
                    q.push(seq[q_pos] as char);
                    q_pos += 1;
                }
                EditOp::Stop => {}
            }
        }
        (q, g)
    }
    pub fn generate(seqs: &[&[u8]], ws: &[f64], config: &Config) -> POA {
        if seqs.is_empty() {
            panic!("Empty string.")
        }
        let seed = (10432940. * ws.iter().sum::<f64>().floor()) as u64 + seqs.len() as u64;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        let max_len = seqs
            .iter()
            .zip(ws.iter())
            .map(|(xs, &w)| if w > 0.001 { xs.len() } else { 0 })
            .max()
            .unwrap_or_else(|| panic!("Empty string."));
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| (&seqs[idx], ws[idx]))
            .filter(|&(_, w)| w > 0.001)
            .fold(POA::default(), |x, (y, w)| {
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
    pub fn generate_w_param<F>(seqs: &[&[u8]], ws: &[f64], ins: i32, del: i32, score: &F) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seed = (10432940. * ws.iter().sum::<f64>().floor()) as u64 + seqs.len() as u64;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        let max_len = seqs
            .iter()
            .zip(ws.iter())
            .map(|(xs, &w)| if w > 0.001 { xs.len() } else { 0 })
            .max()
            .unwrap_or_else(|| panic!("Empty string."));
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| (&seqs[idx], ws[idx]))
            .filter(|&(_, w)| w > 0.001)
            .fold(POA::default(), |x, (y, w)| {
                if x.nodes.len() > 3 * max_len / 2 {
                    x.add_w_param(y, w, ins, del, score).remove_node()
                } else {
                    x.add_w_param(y, w, ins, del, score)
                }
            })
            .remove_node()
            .clean_up()
            .finalize()
    }
    pub fn generate_w_param_simd<F>(
        seqs: &[&[u8]],
        ws: &[f64],
        ins: i32,
        del: i32,
        score: &F,
    ) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seed = (10432940. * ws.iter().sum::<f64>().floor()) as u64 + seqs.len() as u64;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
        let max_len = seqs
            .iter()
            .zip(ws.iter())
            .map(|(xs, &w)| if w > 0.001 { xs.len() } else { 0 })
            .max()
            .unwrap_or_else(|| panic!("Empty string."));
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| (&seqs[idx], ws[idx]))
            .filter(|&(_, w)| w > 0.001)
            .fold(POA::default(), |x, (y, w)| {
                if x.nodes.len() > 3 * max_len / 2 {
                    x.add_w_param_simd(y, w, ins, del, score).remove_node()
                } else {
                    x.add_w_param_simd(y, w, ins, del, score)
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
        let dp = vec![vec![]];
        let route_weight = vec![vec![]];
        let last_elements = vec![];
        Self {
            weight,
            nodes,
            dp,
            route_weight,
            last_elements,
        }
    }
    fn initialize(&mut self, row: usize, column: usize, column_s: i32) {
        let rowvec = std::iter::repeat(MIN).take(column);
        self.dp.iter_mut().for_each(|e| {
            e.clear();
            e.extend(rowvec.clone())
        });
        if self.dp.len() < row {
            for _ in 0..row - self.dp.len() {
                self.dp.push(rowvec.clone().collect());
            }
        }
        for j in 0..column {
            self.dp[0][j] = column_s * j as i32;
        }
        self.route_weight.iter_mut().for_each(|e| e.clear());
        if self.route_weight.len() < row {
            for _ in 0..row - self.route_weight.len() {
                self.route_weight.push(Vec::with_capacity(column));
            }
        }
        let rowvec = std::iter::repeat(0.).take(column);
        for row_vec in self.route_weight.iter_mut() {
            row_vec.extend(rowvec.clone());
        }
        self.last_elements.clear();
        self.last_elements.push(MIN);
    }
    pub fn align_simd<F>(&mut self, seq: &[u8], ins: i32, del: i32, score: F) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        use packed_simd::f32x8 as f32s;
        use packed_simd::i32x8 as i32s;
        const LANE: usize = 8;
        let profile: Vec<Vec<_>> = b"ACGT"
            .iter()
            .map(|&base| seq.iter().map(|&b| score(base, b)).collect())
            .collect();
        let edges: Vec<Vec<_>> = self
            .reverse_edges()
            .iter()
            .map(|edges| {
                if !edges.is_empty() {
                    edges.iter().map(|&p| p + 1).collect()
                } else {
                    vec![0]
                }
            })
            .collect();
        let deletions = i32s::splat(del);
        // -----> query position ---->
        // 0 8 8 8 88 8 8 8 88
        // 0
        // 0
        // |
        // |
        // Graph position
        // |
        // v
        let (row, column) = (self.nodes.len() + 1, seq.len() + 1);
        // Initialazation.
        // self.initialize(row, column, ins);
        // let mut dp = &mut self.dp;
        // let mut last_elements = &mut self.last_elements;
        // let mut route_weight = &mut self.route_weight;
        let mut dp = vec![vec![MIN; column]; row];
        // The last elements in each row. Used at the beggining of the trace-back.
        let mut last_elements = vec![MIN];
        // Weight to break tie. [i][j] is the total weight to (i,j) element.
        let mut route_weight = vec![vec![0.; column]; row];
        for j in 0..column {
            dp[0][j] = ins * j as i32;
            route_weight[0][j] = j as f32;
        }
        for i in 0..row {
            dp[i][0] = 0;
        }
        // The leftmost element used in the match transition.
        for i in 1..row {
            let bs = base_table::BASE_TABLE[self.nodes[i - 1].base() as usize];
            let nw = self.nodes[i - 1].weight() as f32;
            let node_weight = f32s::splat(nw);
            for &p in &edges[i - 1] {
                // Update by SIMD instructions.
                for j in 0..seq.len() / LANE {
                    let (start, end) = (LANE * j + 1, LANE * (j + 1) + 1);
                    // Location to be updated.
                    let current = i32s::from_slice_unaligned(&dp[i][start..end]);
                    let current_weight = f32s::from_slice_unaligned(&route_weight[i][start..end]);
                    // Update for deletion state.
                    let deletion = i32s::from_slice_unaligned(&dp[p][start..end]) + deletions;
                    let deletion_weight = f32s::from_slice_unaligned(&route_weight[p][start..end]);
                    let mask = deletion.gt(current)
                        | (deletion.eq(current) & deletion_weight.gt(current_weight));
                    let current = deletion.max(current);
                    let current_weight = mask.select(deletion_weight, current_weight);
                    // Update for match state.
                    let match_s = i32s::from_slice_unaligned(&profile[bs][start - 1..end - 1])
                        + i32s::from_slice_unaligned(&dp[p][start - 1..end - 1]);
                    let match_weight = node_weight
                        + f32s::from_slice_unaligned(&route_weight[p][start - 1..end - 1]);
                    let mask = match_s.gt(current)
                        | (match_s.eq(current) & match_weight.gt(current_weight));
                    let current = match_s.max(current);
                    let current_weight = mask.select(match_weight, current_weight);
                    current.write_to_slice_unaligned(&mut dp[i][start..end]);
                    current_weight.write_to_slice_unaligned(&mut route_weight[i][start..end]);
                }
                // Update by usual updates.
                for j in (seq.len() / LANE) * LANE..seq.len() {
                    let pos = j + 1;
                    let (del, del_weight) = (dp[p][pos] + del, route_weight[p][pos]);
                    let (current, current_weight) = (dp[i][pos], route_weight[i][pos]);
                    if del > current || (del == current && del_weight > current_weight) {
                        dp[i][pos] = del;
                        route_weight[i][pos] = del_weight;
                    }
                    let mat = dp[p][j] + profile[bs][j];
                    let mat_weight = route_weight[p][j] + nw;
                    let (current, current_weight) = (dp[i][pos], route_weight[i][pos]);
                    if mat > current || (mat == current && mat_weight > current_weight) {
                        dp[i][pos] = mat;
                        route_weight[i][pos] = mat_weight;
                    }
                }
            }
            // Insertions would be updated by usual updates due to dependencies.
            for j in 1..column {
                let (ins, ins_weight) = (dp[i][j - 1] + ins, route_weight[i][j - 1] + 1.);
                let (current, current_weight) = (dp[i][j], route_weight[i][j]);
                if ins > current || (ins == current && ins_weight > current_weight) {
                    dp[i][j] = ins;
                    route_weight[i][j] = ins_weight;
                }
            }
            last_elements.push(dp[i][column - 1]);
        }
        // Traceback.
        let mut q_pos = seq.len();
        let (mut g_pos, &score) = last_elements
            .iter()
            .enumerate()
            .max_by(|(a_i, a), (b_i, b)| match a.partial_cmp(&b) {
                Some(std::cmp::Ordering::Equal) => a_i.cmp(&b_i),
                Some(x) => x,
                None => panic!("{}", line!()),
            })
            .unwrap_or_else(|| panic!("{}", line!()));
        let mut operations = vec![];
        'outer: while q_pos > 0 && g_pos > 0 {
            // Determine where dp[g_pos][q_pos] comes from.
            let w = self.nodes[g_pos - 1].weight() as f32;
            let score = dp[g_pos][q_pos];
            let weight = route_weight[g_pos][q_pos];
            // Deletion.
            for &p in &edges[g_pos - 1] {
                let (del, del_w) = (dp[p][q_pos] + del, route_weight[p][q_pos]);
                if del == score && del_w == weight {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = p;
                    continue 'outer;
                }
            }
            // Insertion
            let ins = dp[g_pos][q_pos - 1] + ins;
            let ins_w = route_weight[g_pos][q_pos - 1] + 1.;
            if ins == score && ins_w == weight {
                q_pos -= 1;
                operations.push(EditOp::Insertion(0));
                continue 'outer;
            }
            // Match/Mismatch
            let bs = base_table::BASE_TABLE[self.nodes[g_pos - 1].base() as usize];
            for &p in &edges[g_pos - 1] {
                let mat = dp[p][q_pos - 1] + profile[bs][q_pos - 1];
                let mat_w = route_weight[p][q_pos - 1] + w;
                if mat == score && mat_w == weight {
                    operations.push(EditOp::Match(g_pos - 1));
                    g_pos = p;
                    q_pos -= 1;
                    continue 'outer;
                }
            }
            panic!("error. none of choices match the current trace table.");
        }
        while q_pos > 0 {
            operations.push(EditOp::Insertion(0));
            q_pos -= 1;
        }
        operations.reverse();
        (score, operations)
    }
    //  ----> Graph position ----->
    // 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  |
    //  |
    //Qeury position
    //  |
    //  v
    // Semi-Global alignment containing query sequence.
    pub fn align<F>(&mut self, seq: &[u8], ins: i32, del: i32, score: F) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        let mut traceback = Vec::with_capacity(seq.len() + 1);
        let edges = self.reverse_edges();
        let mut prev = vec![0; self.nodes.len() + 1];
        let mut p_weight = vec![0.; self.nodes.len() + 1];
        let tb = vec![EditOp::Stop; self.nodes.len() + 1];
        traceback.push(tb);
        let mut updated = vec![];
        let mut updated_weight = vec![];
        let min = -100_000;
        use std::cmp::Ordering::*;
        for (i, &b) in seq.iter().enumerate() {
            updated.push(ins * (i + 1) as i32);
            updated_weight.push((i + 1) as f64);
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
                        .fold((0, min, 0, min), |(m_ax, m_x, d_ax, d_x), (i, m, d)| {
                            let (m_ax, m_x) = match m_x.cmp(&m) {
                                Greater => (m_ax, m_x),
                                Less => (i, m),
                                Equal if p_weight[m_ax] >= p_weight[i] => (m_ax, m_x),
                                Equal => (i, m),
                            };
                            let (d_ax, d_x) = match d_x.cmp(&d) {
                                Greater => (d_ax, d_x),
                                Less => (i, d),
                                Equal if p_weight[m_ax] >= p_weight[i] => (d_ax, d_x),
                                Equal => (i, d),
                            };
                            (m_ax, m_x, d_ax, d_x)
                        })
                };
                let (i_max, i_weight) = (prev[j + 1] + ins, p_weight[j + 1] + 1.);
                let d_weight = updated_weight[d_argmax];
                let m_weight = p_weight[m_argmax] + self.nodes[j].weight();
                use EditOp::*;
                let (max_weight, argmax, max_score) = match d_max.cmp(&i_max) {
                    Greater => (d_weight, Deletion(d_argmax), d_max),
                    Less => (i_weight, Insertion(0), i_max),
                    Equal if d_weight >= i_weight => (d_weight, Deletion(d_argmax), d_max),
                    Equal => (i_weight, Insertion(0), i_max),
                };
                let (max_weight, argmax, max_score) = match max_score.cmp(&m_max) {
                    Greater => (max_weight, argmax, max_score),
                    Less => (m_weight, Match(m_argmax), m_max),
                    Equal if max_weight >= m_weight => (max_weight, argmax, max_score),
                    Equal => (m_weight, Match(m_argmax), m_max),
                };
                assert_eq!(max_score, d_max.max(m_max).max(i_max));
                updated.push(max_score);
                updated_weight.push(max_weight);
                tb.push(argmax)
            }
            traceback.push(tb);
            std::mem::swap(&mut prev, &mut updated);
            std::mem::swap(&mut p_weight, &mut updated_weight);
            updated.clear();
            updated_weight.clear();
        }
        // Traceback
        // g_pos = position on the graph, q_pos = position on the query
        let mut q_pos = seq.len();
        let mut operations = vec![];
        let (mut g_pos, poa_score) = prev
            .into_iter()
            .enumerate()
            .max_by(|(a_i, a), (b_i, b)| match a.partial_cmp(&b) {
                Some(std::cmp::Ordering::Equal) => a_i.cmp(&b_i),
                Some(x) => x,
                None => panic!("{}", line!()),
            })
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
                    operations.push(EditOp::Insertion(0));
                    q_pos -= 1;
                }
                EditOp::Stop => break,
            }
        }
        operations.reverse();
        (poa_score, operations)
    }
    pub fn add_with(self, seq: &[u8], w: f64, c: &Config) -> Self {
        let ins = (c.p_ins.ln() * 3.).floor() as i32;
        let del = (c.p_del.ln() * 3.).floor() as i32;
        let mat = (-10. * c.p_match.ln() * 3.).floor() as i32;
        let mism = (c.mismatch.ln() * 3.).floor() as i32;
        self.add_w_param(seq, w, ins, del, |x, y| if x == y { mat } else { mism })
    }
    pub fn add(self, seq: &[u8], w: f64) -> Self {
        self.add_w_param(seq, w, -2, -2, |x, y| if x == y { 1 } else { -1 })
    }
    pub fn add_w_param<F>(mut self, seq: &[u8], w: f64, ins: i32, del: i32, score: F) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, w);
        }
        // Alignment
        let (_, traceback) = self.align(seq, ins, del, score);
        self.integrate_alignment(seq, w, traceback)
    }
    pub fn add_w_param_simd<F>(mut self, seq: &[u8], w: f64, ins: i32, del: i32, score: F) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, w);
        }
        // Alignment
        let (_, traceback) = self.align_simd(seq, ins, del, &score);
        self.integrate_alignment(seq, w, traceback)
    }

    fn integrate_alignment(mut self, seq: &[u8], w: f64, traceback: TraceBack) -> Self {
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
                    if !edges[from].contains(&to) {
                        edges[from].push(to);
                    }
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
                    if !edges[to].contains(&from) {
                        edges[to].push(from);
                    }
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
