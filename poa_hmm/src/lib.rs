#![feature(test)]
extern crate log;
extern crate packed_simd;
extern crate rand;
extern crate rand_xoshiro;
#[cfg(test)]
extern crate rayon;
#[cfg(test)]
extern crate test;
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
// const LAMBDA_INS: f64 = 0.04;
// const LAMBDA_MATCH: f64 = 0.04;
const LAMBDA_INS: f64 = 0.09;
const LAMBDA_MATCH: f64 = 0.09;
const THR: f64 = 0.4;
pub const DEFAULT_LK: f64 = -150.;
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
    dfs_flag: Vec<usize>,
    dfs_stack: Vec<usize>,
    mapping: Vec<usize>,
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
    pub fn weight(&self) -> f64 {
        self.weight
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
    pub fn generate_by(seqs: &[&[u8]], ws: &[f64], c: &Config) -> POA {
        let ins = (c.p_ins.ln() * 3.).floor() as i32;
        let del = (c.p_del.ln() * 3.).floor() as i32;
        let mat = (-10. * c.p_match.ln() * 3.).floor() as i32;
        let mism = (c.mismatch.ln() * 3.).floor() as i32;
        let score = |x, y| if x == y { mat } else { mism };
        Self::default().update_auto(seqs, ws, (ins, del, &score))
    }
    pub fn from_vec<F>(seqs: &[Vec<u8>], ws: &[f64], parameters: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
        Self::generate(&seqs, ws, parameters)
    }
    pub fn generate<F>(seqs: &[&[u8]], ws: &[f64], parameters: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seed = 99_999_111 * ((ws.iter().sum::<f64>().floor()) as u64);
        Self::default().update(seqs, ws, parameters, seed)
    }
    pub fn update_auto<F>(self, seqs: &[&[u8]], ws: &[f64], params: (i32, i32, &F)) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let seed = 99_999_111 * ((ws.iter().sum::<f64>().floor()) as u64);
        self.update(seqs, ws, params, seed)
    }
    pub fn update<F>(self, seqs: &[&[u8]], ws: &[f64], params: (i32, i32, &F), s: u64) -> POA
    where
        F: Fn(u8, u8) -> i32,
    {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(s);
        if seqs.is_empty() || ws.iter().all(|&w| w <= 0.001) {
            return self;
        }
        let max_len = seqs
            .iter()
            .zip(ws)
            .filter(|&(_, &w)| w > 0.001)
            .map(|(s, _)| s.len())
            .max()
            .unwrap_or(0);
        rand::seq::index::sample(&mut rng, seqs.len(), seqs.len())
            .into_iter()
            .map(|idx| (&seqs[idx], ws[idx]))
            .filter(|&(_, w)| w > 0.001)
            .fold(self, |x, (y, w)| {
                if x.nodes.len() > 3 * max_len / 2 {
                    x.add(y, w, params).remove_node(THR)
                } else {
                    x.add(y, w, params)
                }
            })
            .remove_node(THR)
            .finalize()
    }
    pub fn generate_uniform(seqs: &[&[u8]]) -> POA {
        let ws = vec![1.; seqs.len()];
        POA::generate_by(seqs, &ws, &DEFAULT_CONFIG)
    }
    pub fn generate_vec(seqs: &[Vec<u8>]) -> POA {
        let ws = vec![1.; seqs.len()];
        let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
        POA::generate_by(&seqs, &ws, &DEFAULT_CONFIG)
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
        let (dfs_flag, dfs_stack, mapping) = (vec![], vec![], vec![]);
        Self {
            weight,
            nodes,
            dfs_stack,
            dfs_flag,
            mapping,
        }
    }
    pub fn align<F>(&self, seq: &[u8], param: (i32, i32, F)) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        let (mut dp, mut profile) = (vec![], vec![]);
        self.align_with_buf(seq, param, (&mut dp, &mut profile))
    }
    pub fn align_with_buf<F>(
        &self,
        seq: &[u8],
        (ins, del, score): (i32, i32, F),
        (dp, profile): (&mut Vec<i32>, &mut Vec<i32>),
    ) -> (i32, TraceBack)
    where
        F: Fn(u8, u8) -> i32,
    {
        // -----> query position ---->
        // 0 8 8 8 88 8 8 8 88
        // 0
        // 0
        // |
        // |
        // Graph position
        // |
        // v
        use packed_simd::f32x8 as f32s;
        use packed_simd::i32x8 as i32s;
        const LANE: usize = 8;
        let row = self.nodes.len() + 1;
        let column = if seq.len() % LANE == 0 {
            seq.len() + 1
        } else {
            seq.len() + (LANE - (seq.len() % LANE)) + 1
        };
        profile.clear();
        dp.clear();
        assert!((column - 1) % LANE == 0);
        for &base in b"ACGT" {
            for i in 0..column {
                if i < seq.len() {
                    profile.push(score(base, seq[i]))
                } else {
                    profile.push(0)
                }
            }
        }
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
        // Initialazation.
        for _ in 0..column * row {
            dp.push(std::i32::MIN);
        }
        // The last elements in each row. Used at the beggining of the trace-back.
        let mut last_elements = vec![std::i32::MIN];
        // Weight to break tie. [i*column + j] is the total weight to (i,j) element.
        let mut route_weight = vec![0.; column * row];
        for j in 0..column {
            dp[j] = ins * j as i32;
            route_weight[j] = j as f32;
        }
        for i in 0..row {
            dp[i * column] = 0;
        }
        let deletions = i32s::splat(del);
        for i in 1..row {
            let bs = base_table::BASE_TABLE[self.nodes[i - 1].base() as usize];
            let nw = self.nodes[i - 1].weight() as f32;
            let node_weight = f32s::splat(nw);
            // Update by SIMD instructions.
            for &p in &edges[i - 1] {
                for lj in (0..column - 1).step_by(LANE) {
                    let start = i * column + lj + 1;
                    let end = start + LANE;
                    let prev_start = p * column + lj + 1;
                    let prev_end = prev_start + LANE;
                    let current = i32s::from_slice_unaligned(&dp[start..end]);
                    let current_weight = f32s::from_slice_unaligned(&route_weight[start..end]);
                    // Update for deletion state.
                    let del_score =
                        i32s::from_slice_unaligned(&dp[prev_start..prev_end]) + deletions;
                    let del_weight =
                        f32s::from_slice_unaligned(&route_weight[prev_start..prev_end]);
                    let mask = current.gt(del_score)
                        | (current.eq(del_score) & current_weight.ge(del_weight));
                    let current = mask.select(current, del_score);
                    let current_weight = mask.select(current_weight, del_weight);
                    // Update for match state.
                    let (seq_s, seq_e) = (column * bs + lj, column * bs + lj + LANE);
                    let mat_score = i32s::from_slice_unaligned(&profile[seq_s..seq_e])
                        + i32s::from_slice_unaligned(&dp[prev_start - 1..prev_end - 1]);
                    let mat_weight = node_weight
                        + f32s::from_slice_unaligned(&route_weight[prev_start - 1..prev_end - 1]);
                    let mask = current.gt(mat_score)
                        | (current.eq(mat_score) & current_weight.gt(mat_weight));
                    mask.select(current, mat_score)
                        .write_to_slice_unaligned(&mut dp[start..end]);
                    mask.select(current_weight, mat_weight)
                        .write_to_slice_unaligned(&mut route_weight[start..end]);
                }
            }
            // Insertions would be updated by usual updates due to dependencies.
            for j in 0..seq.len() {
                let pos = i * column + j + 1;
                let current = dp[pos];
                let current_weight = route_weight[pos];
                let ins = dp[pos - 1] + ins;
                let ins_weight = route_weight[pos - 1] + 1.;
                if ins > current || (ins == current && ins_weight > current_weight) {
                    dp[pos] = ins;
                    route_weight[pos] = ins_weight;
                }
            }
            last_elements.push(dp[i * column + seq.len()]);
        }
        // Traceback.
        let mut q_pos = seq.len();
        let (mut g_pos, &score) = last_elements
            .iter()
            .enumerate()
            .max_by(|(a_i, a), (b_i, b)| match a.partial_cmp(&b).unwrap() {
                std::cmp::Ordering::Equal => a_i.cmp(&b_i),
                x => x,
            })
            .unwrap_or_else(|| panic!("{}", line!()));
        assert_eq!(dp[g_pos * column + q_pos], score);
        let mut operations = vec![];
        'outer: while q_pos > 0 && g_pos > 0 {
            let w = self.nodes[g_pos - 1].weight() as f32;
            let score = dp[g_pos * column + q_pos];
            let weight = route_weight[g_pos * column + q_pos];
            // Deletion.
            for &p in &edges[g_pos - 1] {
                let pos = p * column + q_pos;
                let (del, del_w) = (dp[pos] + del, route_weight[pos]);
                if del == score && ((del_w - weight) / del_w.max(weight)).abs() < 0.000_1 {
                    operations.push(EditOp::Deletion(g_pos - 1));
                    g_pos = p;
                    continue 'outer;
                }
            }
            // Insertion
            let ins = dp[g_pos * column + q_pos - 1] + ins;
            let ins_w = route_weight[g_pos * column + q_pos - 1] + 1.;
            if ins == score && ((ins_w - weight) / ins_w.max(weight)).abs() < 0.000_1 {
                q_pos -= 1;
                operations.push(EditOp::Insertion(0));
                continue 'outer;
            }
            // Match/Mismatch
            let bs = base_table::BASE_TABLE[self.nodes[g_pos - 1].base() as usize];
            for &p in &edges[g_pos - 1] {
                let mat = dp[p * column + q_pos - 1] + profile[bs * column + q_pos - 1];
                let mat_w = route_weight[p * column + q_pos - 1] + w;
                if mat == score && ((mat_w - weight) / mat_w.max(weight)).abs() < 0.000_1 {
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
    pub fn add_default(self, seq: &[u8], w: f64) -> Self {
        self.add(seq, w, (-2, -2, &|x, y| if x == y { 1 } else { -1 }))
    }
    pub fn add<F>(self, seq: &[u8], w: f64, parameters: (i32, i32, &F)) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        let mut buf1 = vec![];
        let mut buf2 = vec![];
        self.add_with_buf(seq, w, parameters, (&mut buf1, &mut buf2))
    }
    pub fn add_with_buf<F>(
        self,
        seq: &[u8],
        w: f64,
        ps: (i32, i32, &F),
        buffers: (&mut Vec<i32>, &mut Vec<i32>),
    ) -> Self
    where
        F: Fn(u8, u8) -> i32,
    {
        if self.weight < SMALL || self.nodes.is_empty() {
            return Self::new(seq, w);
        }
        let (_, traceback) = self.align_with_buf(seq, ps, buffers);
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
        // assert!(self.nodes.iter().all(|node| node.weight() > 0.));
        self.topological_sort()
    }
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
    fn finalize(mut self) -> Self {
        let bases: Vec<_> = self.nodes.iter().map(|e| e.base).collect();
        self.nodes.iter_mut().for_each(|e| e.finalize(&bases));
        self
    }
}
