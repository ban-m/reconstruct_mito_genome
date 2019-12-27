#![feature(test)]
//! A tiny implementation for de Bruijn graph with Hidden Markov model.
//! Currently, this implementation is minimal. In other words, it exposes only one struct with just two methods:
//! [DeBruijnGraphHiddenMarkovModel] -- Yes, it is too long -- ,[constructor](DeBruijnGraphHiddenMarkovModel::new),
//! and the [forwardalgorithm](DeBruijnGraphHiddenMarkovModel::forward)
//! to calculate the probability this graph would generate the given observation.
//! As a shorthand for the vary long name, I also supply [DBGHMM] as a alias for [DeBruijnGraphHiddenMarkovModel].
extern crate edlib_sys;
extern crate env_logger;
#[allow(unused_imports)]
#[macro_use]
extern crate log;
extern crate packed_simd;
extern crate rand;
extern crate rand_xoshiro;
extern crate test;
// Whether or not to use 'pseudo count' in the out-dgree.
const PSEUDO_COUNT: f64 = 1.0;
const THR_ON: bool = true;
// This 2.5 is good for prediction for a new query,
// but very bad for exisiting(i.e., the reads used as model) query.
// const THR: f64 = 2.0;
const THR: f64 = 2.0;
// This is tuned for clustering.
// const THR: f64 = 3.;
// const WEIGHT_THR: f64 = 2.0;
// const LOW_LIKELIHOOD: f64 = -100_000.;
const SCALE: f64 = 3.;
const PRIOR_FACTOR: f64 = 0.05;
mod find_union;
pub mod gen_sample;
mod kmer;
use kmer::Kmer;
// This setting is determined by experimentally.
mod config;
pub use config::*;
use packed_simd::f64x4 as f64s;
pub type DBGHMM = DeBruijnGraphHiddenMarkovModel;
mod factory;
pub use factory::Factory;
mod base_table;

#[derive(Debug, Clone, Default)]
pub struct DeBruijnGraphHiddenMarkovModel {
    nodes: Vec<Kmer>,
    k: usize,
    // roughtly the same as coverage.
    weight: f64,
}

impl std::fmt::Display for DBGHMM {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edges = self
            .nodes
            .iter()
            .map(|kmer| kmer.edges.iter().filter(|e| e.is_some()).count())
            .sum::<usize>();
        write!(
            f,
            "K:{}\tNodes:{}\tEdges:{}\tWeight:{:.3}",
            self.k,
            self.nodes.len(),
            edges,
            self.weight
        )
    }
}

impl DeBruijnGraphHiddenMarkovModel {
    pub fn from(nodes: Vec<Kmer>, k: usize, weight: f64) -> Self {
        Self { nodes, k, weight }
    }
    pub fn new(dataset: &[Vec<u8>], k: usize) -> Self {
        let mut f = Factory::new();
        f.generate(dataset, k)
    }
    pub fn new_from_ref(dataset: &[&[u8]], k: usize) -> Self {
        let mut f = Factory::new();
        f.generate_from_ref(dataset, k)
    }
    pub fn new_with_weight(dataset: &[&[u8]], ws: &[f64], k: usize) -> Self {
        let mut f = Factory::new();
        f.generate_with_weight(dataset, ws, k)
    }
    pub fn new_with_weight_prior(dataset: &[&[u8]], ws: &[f64], k: usize) -> Self {
        let mut f = Factory::new();
        let mut buf = vec![];
        f.generate_with_weight_prior(dataset, ws, k, &mut buf)
    }
    // Calc hat. Return (c,d)
    fn is_non_zero(idx: usize, xs: &[f64]) -> bool {
        xs[idx * 3..(idx + 1) * 3].iter().any(|&e| e > 0.00001)
    }
    fn sum(xs: &[f64]) -> f64 {
        assert!(xs.len() % 4 == 0);
        xs.chunks_exact(f64s::lanes())
            .map(f64s::from_slice_unaligned)
            .sum::<f64s>()
            .sum()
    }
    fn mul(xs: &mut [f64], y: f64) {
        assert!(xs.len() % 4 == 0);
        let ys = f64s::splat(y);
        xs.chunks_exact_mut(f64s::lanes()).for_each(|xs| {
            let packed = f64s::from_slice_unaligned(xs) * ys;
            packed.write_to_slice_unaligned(xs);
        });
    }
    fn update(
        &self,
        updates: &mut [f64],
        prev: &[f64],
        base: u8,
        config: &Config,
        edges: &[Vec<(usize, usize)>],
    ) -> (f64, f64) {
        debug_assert!((1. - prev.iter().sum::<f64>()).abs() < 0.001);
        // Alignemnt:[mat,ins,del, mat,ins,del,, mat,....,del]
        for (idx, from) in self
            .nodes
            .iter()
            .enumerate()
            .filter(|&(idx, _)| Self::is_non_zero(idx, prev))
        {
            let node = 3 * idx;
            let prob = prev[node] * config.p_match
                + prev[node + 1] * (1. - config.p_extend_ins)
                + prev[node + 2] * (1. - config.p_extend_del - config.p_del_to_ins);
            let del_to_ins = prev[node + 2] * config.p_del_to_ins;
            let ins = if !edges[idx].is_empty() {
                prev[node] * config.p_ins + prev[node + 1] * config.p_extend_ins
            } else {
                prev[node] + prev[node + 1] + prev[node + 2]
            };
            updates[node + 1] += ins * from.insertion(base);
            // From -> To with base b
            edges[idx].iter().for_each(|&(i, to)| {
                updates[3 * to] += prob * from.to(i) * self.nodes[to].prob(base, config);
                updates[3 * to + 1] += del_to_ins * from.to(i) * self.nodes[to].insertion(base);
            });
        }
        let d = Self::sum(updates).recip();
        // Updates -> tilde
        Self::mul(updates, d);
        // Update `Del` states.
        for (idx, edges) in edges.iter().enumerate() {
            if Self::is_non_zero(idx, updates) {
                let from = &updates[3 * idx..3 * (idx + 1)];
                let pushing_weight: f64 = from[0] * config.p_del + from[2] * config.p_extend_del;
                for (_, to) in edges.iter() {
                    updates[3 * to + 2] += pushing_weight;
                }
            }
        }
        let c = Self::sum(&updates).recip();
        Self::mul(updates, c);
        (c, d)
    }
    // ln sum_i exp(x[i])
    fn logsumexp(xs: &[f64]) -> f64 {
        let max = xs
            .iter()
            .fold(std::f64::MIN, |x, &y| if x < y { y } else { x });
        max + xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln()
    }
    fn global_alns(&self, tip: &[u8], config: &Config) -> Vec<f64> {
        let (del, ins) = (config.p_del.ln(), config.p_ins.ln());
        let mism = (config.mismatch / 3.).ln();
        let mat = (1. - config.mismatch).ln();
        let mut prev = vec![0.; self.k + 1];
        let mut next = vec![0.; self.k + 1];
        self.nodes
            .iter()
            .map(|e| Self::global_aln(&e.kmer, tip, mat, mism, del, ins, &mut prev, &mut next))
            .collect()
    }
    #[inline]
    fn global_aln(
        xs: &[u8],
        ys: &[u8],
        mat: f64,
        mism: f64,
        del: f64,
        ins: f64,
        prev: &mut Vec<f64>,
        next: &mut Vec<f64>,
    ) -> f64 {
        let small = std::f64::MIN;
        let k = xs.len();
        for i in 0..=k {
            prev[i] = small;
        }
        prev[0] = 0.;
        for i in 0..k {
            next[0] = small;
            for j in 0..k {
                let match_score = prev[j] + if xs[i] == ys[j] { mat } else { mism };
                next[j + 1] = (next[j] + del).max(prev[j + 1] + ins).max(match_score);
            }
            std::mem::swap(prev, next);
        }
        prev[k]
    }
    // Return edit operations.
    #[inline]
    fn global_aln_ops(xs: &[u8], ys: &[u8], mat: f64, mism: f64, del: f64, ins: f64) -> Vec<u8> {
        let mut dp = vec![vec![std::f64::MIN; ys.len() + 1]; xs.len() + 1];
        let mut traceback = vec![vec![0; ys.len() + 1]; xs.len() + 1];
        // DP fill.
        dp[0][0] = 0.;
        for i in 0..xs.len() {
            for j in 0..ys.len() {
                // Fill dp[i+1][j+1]
                let match_score = dp[i][j] + if xs[i] == ys[j] { mat } else { mism };
                let del_score = dp[i][j + 1] + del;
                let ins_score = dp[i + 1][j] + ins;
                if del_score < match_score && ins_score < match_score {
                    dp[i + 1][j + 1] = match_score;
                    traceback[i + 1][j + 1] = 0;
                } else if match_score < ins_score && del_score < ins_score {
                    dp[i + 1][j + 1] = ins_score;
                    traceback[i + 1][j + 1] = 1;
                } else {
                    dp[i + 1][j + 1] = del_score;
                    traceback[i + 1][j + 1] = 2;
                }
            }
        }
        // Traceback.
        let mut ops = vec![];
        let (mut row, mut column) = (xs.len(), ys.len());
        while row > 0 && column > 0 {
            let op = traceback[row][column];
            ops.push(op);
            if op == 0 || op == 1 {
                column -= 1;
            }
            if op == 0 || op == 2 {
                row -= 1;
            }
        }
        ops.reverse();
        ops
    }
    fn initialize(&self, tip: &[u8], config: &Config) -> (f64, f64, Vec<f64>) {
        let aln_scores = self.global_alns(tip, config);
        let sum = Self::logsumexp(&aln_scores);
        let mut hat = Vec::with_capacity(aln_scores.len() * 3);
        for init in aln_scores.into_iter().map(|x| (x - sum).exp()) {
            hat.push(init);
            hat.push(0.);
            hat.push(0.);
        }
        hat.extend(std::iter::repeat(0.).take(4 - hat.len() % 4));
        assert!(hat.len() % 4 == 0);
        let (minus_ln_c, minus_ln_d) = (0., 0.);
        (minus_ln_c, minus_ln_d, hat)
    }
    /// Return the total weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    #[cfg(target_feature = "sse")]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        // Alignemnts: [mat, ins, del, mat, ins, del, ....]
        let (mut cs, mut ds, mut prev) = self.initialize(&obs[..self.k], config);
        assert!(prev.len() % 4 == 0);
        let mut updated = vec![0.; prev.len()];
        let edges: Vec<Vec<(usize, usize)>> = self
            .nodes
            .iter()
            .map(|e| {
                e.edges
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, e)| e.map(|to| (idx, to)))
                    .collect()
            })
            .collect();
        for &base in &obs[self.k..] {
            updated.iter_mut().for_each(|e| *e = 0.);
            let (c, d) = self.update(&mut updated, &prev, base, config, &edges);
            cs -= c.ln();
            ds -= d.ln();
            std::mem::swap(&mut prev, &mut updated);
        }
        assert!(cs + ds < 0.);
        cs + ds
    }
    fn update_exp(
        &self,
        updates: &mut [f64],
        prev: &[f64],
        base: u8,
        config: &Config,
        edges: &[Vec<usize>],
    ) -> (f64, f64) {
        debug_assert!((1. - prev.iter().sum::<f64>()).abs() < 0.001);
        // Alignemnt:[mat,ins,del, mat,ins,del,, mat,....,del]
        for (dist_idx, dist) in self.nodes.iter().enumerate() {
            let edge_idx = base_table::BASE_TABLE[dist.last() as usize];
            let node = 3 * dist_idx;
            updates[node] = edges[dist_idx]
                .iter()
                .map(|&src| {
                    let src_node = 3 * src;
                    let prob = prev[src_node] * config.p_match
                        + prev[src_node + 1] * (1. - config.p_extend_ins)
                        + prev[src_node + 2] * (1. - config.p_extend_del - config.p_del_to_ins);
                    prob * self.nodes[src].to(edge_idx)
                })
                .sum::<f64>()
                * dist.prob(base, config);
            updates[node + 1] = if dist.edges.iter().any(|e| e.is_some()) {
                prev[node] * config.p_ins + prev[node + 1] * config.p_extend_ins
            } else {
                prev[node..node + 3].iter().sum::<f64>()
            };
            updates[node + 1] += edges[dist_idx]
                .iter()
                .map(|&src| {
                    let trans = prev[3 * src + 2] * config.p_del_to_ins;
                    trans * self.nodes[src].to(edge_idx)
                })
                .sum::<f64>();
            updates[node + 1] *= dist.insertion(base);
        }
        let d = Self::sum(updates).recip();
        Self::mul(updates, d);
        for (dist_idx, src_nodes) in edges.iter().enumerate() {
            let edge_idx = base_table::BASE_TABLE[self.nodes[dist_idx].last() as usize];
            updates[3 * dist_idx + 2] += src_nodes
                .iter()
                .map(|&src| {
                    let trans = updates[3 * src] * config.p_del
                        + updates[3 * src + 2] * config.p_extend_del;
                    trans * self.nodes[src].to(edge_idx)
                })
                .sum::<f64>();
        }
        let c = Self::sum(&updates).recip();
        Self::mul(updates, c);
        (c, d)
    }
    #[cfg(target_feature = "sse")]
    pub fn forward_exp(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        // Alignemnts: [mat, ins, del, mat, ins, del, ....]
        let (mut cs, mut ds, mut prev) = self.initialize(&obs[..self.k], config);
        assert!(prev.len() % 4 == 0);
        let mut updated = vec![0.; prev.len()];
        let edges = {
            let mut edges = vec![Vec::with_capacity(4); self.nodes.len()];
            for (from, n) in self.nodes.iter().enumerate() {
                for &to in n.edges.iter().filter_map(|e| e.as_ref()) {
                    edges[to].push(from);
                }
            }
            edges
        };
        for &base in &obs[self.k..] {
            updated.iter_mut().for_each(|e| *e = 0.);
            let (c, d) = self.update_exp(&mut updated, &prev, base, config, &edges);
            cs -= c.ln();
            ds -= d.ln();
            std::mem::swap(&mut prev, &mut updated);
        }
        cs + ds
    }
    pub fn node_num(&self) -> usize {
        self.nodes.len()
    }
    pub fn edge_num(&self) -> usize {
        self.nodes
            .iter()
            .map(|kmer| kmer.edges.iter().filter(|e| e.is_some()).count())
            .sum::<usize>()
    }
    pub fn has(&self, xs: &[u8]) -> bool {
        self.nodes.iter().any(|kmer| kmer.kmer == xs)
    }
    pub fn is_connected_fu(&self) -> bool {
        let mut fu = find_union::FindUnion::new(self.nodes.len());
        for (idx, n) in self.nodes.iter().enumerate() {
            for to in n.edges.iter().filter_map(|e| e.as_ref()) {
                fu.unite(idx, *to).unwrap();
            }
        }
        let node: std::collections::HashSet<_> =
            (0..self.nodes.len()).map(|e| fu.find(e).unwrap()).collect();
        node.len() == 1
    }
    pub fn is_connected(&self) -> bool {
        let mut nodes: Vec<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|&(_, e)| e.is_head)
            .collect();
        let mut is_reached: Vec<_> = vec![false; self.nodes.len()];
        'dfs: while !nodes.is_empty() {
            let &(idx, node) = nodes.last().unwrap();
            is_reached[idx] = true;
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                if !is_reached[to] {
                    nodes.push((to, &self.nodes[to]));
                    continue 'dfs;
                }
            }
            nodes.pop().unwrap();
        }
        is_reached.into_iter().all(|e| e)
    }
    // Compute a row of Viterbi algorithm. Return the corresponding traceback row.
    fn viterbi_row(
        &self,
        base: u8,
        prev: &[f64],
        updates: &mut [f64],
        config: &Config,
        edges: &[Vec<usize>],
    ) -> Vec<usize> {
        use std::cmp::Ordering::Equal;
        // Compute the match/insertion state.
        let mut arg_updates = vec![0; updates.len()];
        for (dist_idx, dist) in self.nodes.iter().enumerate() {
            let edge_idx = base_table::BASE_TABLE[dist.last() as usize];
            let node = 3 * dist_idx;
            let emit = dist.prob(base, config).ln();
            let (argmax, max) = edges[dist_idx]
                .iter()
                .map(|&src| {
                    let src_node = 3 * src;
                    let (mat, ins, del) = (
                        prev[src_node] + config.p_match.ln(),
                        prev[src_node + 1] + (1. - config.p_extend_ins).ln(),
                        prev[src_node + 2] + (1. - config.p_extend_del - config.p_del_to_ins).ln(),
                    );
                    let trans = self.nodes[src].to(edge_idx).ln();
                    if del < mat && ins < mat {
                        (src_node, mat + trans + emit)
                    } else if mat < ins && del < ins {
                        (src_node + 1, ins + trans + emit)
                    } else {
                        (src_node + 2, del + trans + emit)
                    }
                })
                .max_by(|e, f| (e.1).partial_cmp(&(f.1)).unwrap_or(Equal))
                .unwrap_or((0, -1000000.));
            updates[node] = max;
            arg_updates[node] = argmax;
            //Ins state
            let (argmax, max) = if dist.edges.iter().any(|e| e.is_some()) {
                let match_to_ins = prev[node] + config.p_ins.ln();
                let ins_extend = prev[node + 1] + config.p_extend_ins.ln();
                if ins_extend < match_to_ins {
                    (node, match_to_ins)
                } else {
                    (node + 1, ins_extend)
                }
            } else {
                let (mat, ins, del) = (prev[node], prev[node + 1], prev[node + 2]);
                if ins < mat && del < mat {
                    (node, mat)
                } else if mat < ins && del < ins {
                    (node + 1, ins)
                } else {
                    (node + 2, del)
                }
            };
            let (argmax_del, max_del) = edges[dist_idx]
                .iter()
                .map(|&src| {
                    let trans = prev[3 * src + 2] + config.p_del_to_ins.ln();
                    (3 * src + 2, trans + self.nodes[src].to(edge_idx).ln())
                })
                .max_by(|e, f| (e.1).partial_cmp(&(f.1)).unwrap_or(Equal))
                .unwrap_or((0, -1000000.));
            if max_del < max {
                updates[node + 1] = max + dist.insertion(base).ln();
                arg_updates[node + 1] = argmax;
            } else {
                updates[node + 1] = max_del + dist.insertion(base).ln();
                arg_updates[node + 1] = argmax_del;
            }
        }
        // Deletion state.
        for (dist_idx, src_nodes) in edges.iter().enumerate() {
            let edge_idx = base_table::BASE_TABLE[self.nodes[dist_idx].last() as usize];
            let (argmax, max) = src_nodes
                .iter()
                .map(|&src| {
                    let mat_del = updates[3 * src] + config.p_del.ln();
                    let del_ext = updates[3 * src + 2] + config.p_extend_del.ln();
                    if del_ext < mat_del {
                        (3 * src, mat_del + self.nodes[src].to(edge_idx).ln())
                    } else {
                        (3 * src + 2, del_ext + self.nodes[src].to(edge_idx).ln())
                    }
                })
                .max_by(|e, f| (e.1).partial_cmp(&(f.1)).unwrap_or(Equal))
                .unwrap_or((0, -1000000.));
            updates[dist_idx * 3 + 2] = max;
            arg_updates[dist_idx * 3 + 2] = argmax;
        }
        arg_updates
    }
    // Viterbi algorithm. Return Vec<(base,state)>.
    pub fn viterbi(&self, obs: &[u8], config: &Config) -> (f64, Vec<(u8, u8)>) {
        assert!(obs.len() > self.k);
        // Alignemnts: [mat, ins, del, mat, ins, del, ....]
        let (_, _, prev) = self.initialize(&obs[..self.k], config);
        let mut prev: Vec<_> = prev
            .into_iter()
            .map(|e| if e == 0. { -10000000. } else { e.ln() })
            .collect();
        let len = prev.len();
        let mut traceback = vec![vec![0; len]];
        let edges = {
            let mut edges = vec![Vec::with_capacity(4); self.nodes.len()];
            for (from, n) in self.nodes.iter().enumerate() {
                for &to in n.edges.iter().filter_map(|e| e.as_ref()) {
                    edges[to].push(from);
                }
            }
            edges
        };
        let mut updates = vec![-100000.; len];
        for &base in &obs[self.k..] {
            let tb = self.viterbi_row(base, &prev, &mut updates, config, &edges);
            traceback.push(tb);
            std::mem::swap(&mut prev, &mut updates);
            updates.iter_mut().for_each(|e| *e = -1000000.);
        }
        // Traceback. `prev` is the last row which was filled.
        // Traceback
        use std::cmp::Ordering::Equal;
        let (mut arg, max) = prev
            .iter()
            .enumerate()
            .max_by(|e, f| (e.1).partial_cmp(&(f.1)).unwrap_or(Equal))
            .unwrap();
        let mut row = traceback.len() - 1;
        let mut result = vec![];
        while row > 0 {
            result.push((self.nodes[arg / 3].last(), (arg % 3) as u8));
            let prev = traceback[row][arg];
            // If deletion, the previous state comes from the same row.
            if prev % 3 != 2 {
                row -= 1;
            }
            arg = prev;
        }

        result.reverse();
        // Tip alignments.
        let (mism, del, ins) = (
            (config.mismatch / 3.).ln(),
            config.p_del.ln(),
            config.p_ins.ln(),
        );
        let mat = (1. - config.mismatch).ln();
        let tip = &self.nodes[arg / 3].kmer;
        let ops = Self::global_aln_ops(tip, &obs[..self.k], mat, mism, del, ins);
        let mut idx = 0;
        let mut base_and_ops: Vec<_> = ops
            .into_iter()
            .map(|op| match op {
                0 | 2 => {
                    idx += 1;
                    (tip[idx - 1], op)
                }
                1 => (b'-', op),
                _ => unreachable!(),
            })
            .collect();
        base_and_ops.append(&mut result);
        (*max, base_and_ops)
    }
    pub fn traverse(&self) -> Vec<u8> {
        use std::cmp::Ordering::Equal;
        let (mut idx, mut node): (usize, &Kmer) = self
            .nodes
            .iter()
            .enumerate()
            .filter(|&(_, ref kmer)| kmer.is_head)
            .max_by(|e, f| ((e.1).tot).partial_cmp(&(f.1).tot).unwrap_or(Equal))
            .unwrap();
        let mut res = node.kmer.clone();
        let mut used = vec![[false; 4]; self.nodes.len()];
        while !node.is_tail {
            let ((next_idx, _), _) = match node.weight[4..]
                .iter()
                .enumerate()
                .zip(used[idx].iter())
                .filter(|&(_, &b)| !b)
                .max_by(|e, f| ((e.0).1).partial_cmp(&(f.0).1).unwrap_or(Equal))
            {
                Some(res) => res,
                None => break,
            };
            used[idx][next_idx] = true;
            idx = match self.nodes[idx].edges[next_idx] {
                Some(res) => res,
                None => break,
            };
            node = &self.nodes[idx];
            res.push(node.last());
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::Xoroshiro128StarStar;
    #[test]
    fn works() {}
    #[test]
    fn initialize() {
        let test = [
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CA".to_vec(),
            b"TTTTTTGTGTGACTGTACGTGACG".to_vec(),
            b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
        ];
        let _ = DBGHMM::new(&test, 12);
    }
    #[test]
    fn forward() {
        let test = [
            b"CAGTGCTAGTCGATGTCA".to_vec(),
            b"CA".to_vec(),
            b"TTTTTTGTGTGACTGTACGTGACG".to_vec(),
            b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
            b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
            b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
        ];
        let mode = DBGHMM::new(&test, 12);
        mode.forward(b"CACACAGCAGTCAGTGCA", &DEFAULT_CONFIG);
    }
    use rand::{seq::SliceRandom, SeedableRng}; // 0.2%, 0.65%, 0.65%.
    struct Profile {
        sub: f64,
        del: f64,
        ins: f64,
    }
    const PROFILE: Profile = Profile {
        sub: 0.03,
        del: 0.05,
        ins: 0.06,
    };
    #[test]
    fn forward_check() {
        //env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template: Vec<_> = (0..30)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..20)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let k = 7;
        let model1 = DBGHMM::new(&model1, k);
        let likelihood1 = model1.forward(&template, &DEFAULT_CONFIG);
        assert!(!likelihood1.is_nan())
    }

    #[test]
    fn random_check() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template: Vec<_> = (0..150)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..50)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..50)
            .map(|_| {
                (0..150)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let k = 7;
        let model1 = DBGHMM::new(&model1, k);
        let model2 = DBGHMM::new(&model2, k);
        let likelihood1 = model1.forward(&template, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&template, &DEFAULT_CONFIG);
        assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
    }
    #[cfg(not(target_feature = "sse"))]
    #[test]
    fn random_check_weight() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template: Vec<_> = (0..150)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..50)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..50)
            .map(|_| {
                (0..150)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let k = 7;
        let weight1: Vec<_> = (0..model1.len())
            .map(|_| 1.)
            .chain((0..model2.len()).map(|_| 0.))
            .collect();
        let weight2: Vec<_> = (0..model1.len())
            .map(|_| 0.)
            .chain((0..model2.len()).map(|_| 1.))
            .collect();
        let dataset: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();
        let mut f = Factory::new();
        let model1 = f.generate_with_weight(&dataset, &weight1, k);
        let model2 = f.generate_with_weight(&dataset, &weight2, k);
        let (a, b, ps) = model1.initialize(&template, &DEFAULT_CONFIG);
        assert!(!a.is_nan() && !b.is_nan());
        assert!(ps.iter().all(|e| !e.fold().is_nan()));
        let mut pps = vec![Node::new(0., 0., 0.); ps.len()];
        let (a, b) = model1.update(&mut pps, &ps, b'A', &DEFAULT_CONFIG);
        assert!(!a.is_nan() && !b.is_nan());
        assert!(pps.into_iter().all(|e| !e.fold().is_nan()));
        let (a, b, ps) = model2.initialize(&template, &DEFAULT_CONFIG);
        assert!(ps.into_iter().all(|e| !e.fold().is_nan()));
        assert!(!a.is_nan() && !b.is_nan());
        let likelihood1 = model1.forward(&template, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&template, &DEFAULT_CONFIG);
        assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
    }
    #[test]
    fn has_true_edge_test() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12_122_432);
        let k = 6;
        for num in 10..200 {
            let template1: Vec<_> = (0..150)
                .filter_map(|_| bases.choose(&mut rng))
                .copied()
                .collect();
            let dataset: Vec<Vec<_>> = (0..num)
                .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
                .collect();
            let weight = vec![1.; num];
            let kmers: std::collections::HashSet<Vec<u8>> = dataset
                .iter()
                .flat_map(|e| e.windows(k))
                .map(|e| e.to_vec())
                .collect();
            let dataset: Vec<_> = dataset.iter().map(|e| e.as_slice()).collect();
            let mut f = Factory::new();
            let model1 = f.generate_with_weight_prior(&dataset, &weight, k, &mut vec![]);
            eprintln!("{}", model1);
            eprintln!("{}", String::from_utf8_lossy(&template1));
            eprintln!("{}", num);
            let num = template1
                .windows(k)
                .filter(|&kmer| !model1.nodes.iter().any(|node| node.kmer == kmer))
                .inspect(|kmer| eprintln!("{}", String::from_utf8_lossy(kmer)))
                .count();
            eprintln!("{}", num);
            if num > 0 {
                for d in dataset {
                    eprintln!("{}", String::from_utf8_lossy(d));
                }
            }
            for kmer in template1.windows(k).filter(|&kmer| kmers.contains(kmer)) {
                assert!(model1.nodes.iter().any(|node| node.kmer == kmer));
            }
        }
    }
    #[test]
    fn hard_test() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template1: Vec<_> = (0..150)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let p = Profile {
            sub: 0.005,
            ins: 0.005,
            del: 0.005,
        };
        let template2 = introduce_randomness(&template1, &mut rng, &p);
        let model1: Vec<Vec<_>> = (0..50)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..50)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let k = 7;
        let model1 = DBGHMM::new(&model1, k);
        let model2 = DBGHMM::new(&model2, k);
        for _i in 0..20 {
            let likelihood1 = model1.forward(&template1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&template1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
        }
        for _i in 0..20 {
            let likelihood1 = model1.forward(&template2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&template2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
        }
    }
    #[test]
    fn mix_test_prior() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template1: Vec<_> = (0..150)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let p = Profile {
            sub: 0.003,
            ins: 0.003,
            del: 0.003,
        };
        let template2 = introduce_randomness(&template1, &mut rng, &p);
        let model1: Vec<Vec<_>> = (0..25)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..25)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let dataset: Vec<_> = model1
            .iter()
            .map(|e| e.as_slice())
            .chain(model2.iter().map(|e| e.as_slice()))
            .collect();
        let mut f = Factory::new();
        let mut buf = vec![];
        let k = 6;
        let weight1 = vec![vec![0.8; 25], vec![0.2; 25]].concat();
        let weight2 = vec![vec![0.2; 25], vec![0.8; 25]].concat();
        let model1 = f.generate_with_weight_prior(&dataset, &weight1, k, &mut buf);
        let model2 = f.generate_with_weight_prior(&dataset, &weight2, k, &mut buf);
        let num = 50;
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template1, &mut rng, &PROFILE);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                lk1 > lk2
            })
            .count();
        // assert!(lk1 > lk2, "{},{},{}", lk1, lk2, i);
        //}
        assert!(correct >= num * 4 / 5, "{}", correct);
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template2, &mut rng, &PROFILE);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                lk1 < lk2
                //assert!(lk1 < lk2, "{},{},{}", lk1, lk2, i);
            })
            .count();
        assert!(correct >= num * 4 / 5, "{}", correct);
    }
    #[test]
    fn single_error_test() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1234565);
        let coverage = 50;
        let start = 7;
        let mut f = Factory::new();
        let results: Vec<_> = (start..coverage)
            .step_by(2)
            .map(|cov| {
                let template1: Vec<_> = (0..150)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 1, 0, 0);
                let sub = check(&template1, &template2, &mut rng, cov, &mut f);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 1, 0);
                let del = check(&template1, &template2, &mut rng, cov, &mut f);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 0, 1);
                let ins = check(&template1, &template2, &mut rng, cov, &mut f);
                (cov, (sub, del, ins))
            })
            .collect();
        let (sub, del, ins) = results
            .iter()
            .fold((0, 0, 0), |(x, y, z), &(_, (a, b, c))| {
                (x + a, y + b, z + c)
            });
        for (cov, res) in results {
            eprintln!("Cov:{},Sub:{},Del:{},Ins:{}", cov, res.0, res.1, res.2);
        }
        eprintln!("Tot:{}", (start..coverage).step_by(2).count() * 100);
        eprintln!("Sub:{},Del:{},Ins:{}", sub, del, ins);
        assert!(false);
    }
    fn check<R: rand::Rng>(
        t1: &[u8],
        t2: &[u8],
        rng: &mut R,
        cov: usize,
        f: &mut Factory,
    ) -> usize {
        let model1: Vec<_> = (0..cov)
            .map(|_| introduce_randomness(&t1, rng, &PROFILE))
            .collect();
        let model2: Vec<_> = (0..cov)
            .map(|_| introduce_randomness(&t2, rng, &PROFILE))
            .collect();
        let seqs: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();
        let weight1 = vec![vec![1.; cov], vec![0.; cov]].concat();
        let weight2 = vec![vec![0.; cov], vec![1.; cov]].concat();
        let k = 6;
        let m1 = f.generate_with_weight(&seqs, &weight1, k);
        let m2 = f.generate_with_weight(&seqs, &weight2, k);
        eprintln!("{}\t{}\t{}", cov, m1, m2);
        let correct = (0..100)
            .filter(|e| {
                if e % 2 == 0 {
                    let q = introduce_randomness(&t1, rng, &PROFILE);
                    // let d1: u32 = model1.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // let d2: u32 = model2.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // d1 < d2
                    m1.forward(&q, &DEFAULT_CONFIG) > m2.forward(&q, &DEFAULT_CONFIG)
                } else {
                    let q = introduce_randomness(&t2, rng, &PROFILE);
                    // let d1: u32 = model1.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // let d2: u32 = model2.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // d1 > d2
                    m1.forward(&q, &DEFAULT_CONFIG) < m2.forward(&q, &DEFAULT_CONFIG)
                }
            })
            .count();
        correct
    }

    #[test]
    fn single_error_prior_test() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1234565);
        let coverage = 80;
        let start = 20;
        let mut f = Factory::new();
        let results: Vec<_> = (start..coverage)
            .step_by(2)
            .map(|cov| {
                let template1: Vec<_> = (0..150)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 1, 0, 0);
                let sub = check_prior(&template1, &template2, &mut rng, cov, &mut f);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 1, 0);
                let del = check_prior(&template1, &template2, &mut rng, cov, &mut f);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 0, 1);
                let ins = check_prior(&template1, &template2, &mut rng, cov, &mut f);
                (cov, (sub, del, ins))
            })
            .collect();
        let (sub, del, ins) = results
            .iter()
            .fold((0, 0, 0), |(x, y, z), &(_, (a, b, c))| {
                (x + a, y + b, z + c)
            });
        for (cov, res) in results {
            eprintln!("Cov:{},Sub:{},Del:{},Ins:{}", cov, res.0, res.1, res.2);
        }
        eprintln!("Sub:{},Del:{},Ins:{}", sub, del, ins);
        eprintln!("Tot:{}", (start..coverage).step_by(2).count() * 100);
        assert!(false);
    }
    fn check_prior<R: rand::Rng>(
        t1: &[u8],
        t2: &[u8],
        rng: &mut R,
        cov: usize,
        f: &mut Factory,
    ) -> usize {
        let mut buf = vec![];
        let model1: Vec<_> = (0..cov)
            .map(|_| introduce_randomness(&t1, rng, &PROFILE))
            .collect();
        let model2: Vec<_> = (0..cov)
            .map(|_| introduce_randomness(&t2, rng, &PROFILE))
            .collect();
        let seqs: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();
        let weight1 = vec![vec![1.; cov], vec![0.; cov]].concat();
        let weight2 = vec![vec![0.; cov], vec![1.; cov]].concat();
        let k = 6;
        let m1 = f.generate_with_weight_prior(&seqs, &weight1, k, &mut buf);
        let m2 = f.generate_with_weight_prior(&seqs, &weight2, k, &mut buf);
        eprintln!("{}\t{}\t{}", cov, m1, m2);
        let correct = (0..100)
            .filter(|e| {
                if e % 2 == 0 {
                    let q = introduce_randomness(&t1, rng, &PROFILE);
                    m1.forward_exp(&q, &DEFAULT_CONFIG) > m2.forward_exp(&q, &DEFAULT_CONFIG)
                // let d1: u32 = model1.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                // let d2: u32 = model2.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                // d1 < d2
                } else {
                    let q = introduce_randomness(&t2, rng, &PROFILE);
                    m1.forward_exp(&q, &DEFAULT_CONFIG) < m2.forward_exp(&q, &DEFAULT_CONFIG)
                    // let d1: u32 = model1.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // let d2: u32 = model2.iter().map(|e| edlib_sys::global_dist(&q, e)).sum();
                    // d1 > d2
                }
            })
            .count();
        correct
    }
    #[test]
    fn low_coverage_test() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(121212);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        for m in &model1 {
            eprintln!("1:{}", String::from_utf8_lossy(m));
        }
        for m in &model2 {
            eprintln!("2:{}", String::from_utf8_lossy(m));
        }
        let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
        let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
        eprintln!("1:{}", String::from_utf8_lossy(&test1));
        eprintln!("2:{}", String::from_utf8_lossy(&test2));
        let model1 = DBGHMM::new(&model1, k);
        eprintln!("Model1:{}", model1);
        let model2 = DBGHMM::new(&model2, k);
        eprintln!("Model2:{}", model2);
        {
            let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
        {
            let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
        // assert!(false);
    }

    #[test]
    fn low_coverage_weighted_test() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(121212);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
        let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
        let dataset: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();

        let weight1 = vec![vec![1.; 10], vec![0.; 10]].concat();
        let weight2 = vec![vec![0.; 10], vec![1.; 10]].concat();
        let mut f = Factory::new();
        let model1 = f.generate_with_weight(&dataset, &weight1, k);
        eprintln!("Model1:{}", model1);
        let model2 = f.generate_with_weight(&dataset, &weight2, k);
        eprintln!("Model2:{}", model2);
        {
            let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
        {
            let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
        // assert!(false);
    }

    #[test]
    fn high_coverage_test() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(121212332);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..200)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..200)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
        let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
        {
            let model1 = DBGHMM::new(&model1, k);
            let model2 = DBGHMM::new(&model2, k);

            let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
            let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
        }
        let dataset: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();
        let weight1 = vec![vec![1.; 200], vec![0.; 200]].concat();
        let weight2 = vec![vec![0.; 200], vec![1.; 200]].concat();
        let mut f = Factory::new();
        let model1 = f.generate_with_weight(&dataset, &weight1, k);
        eprintln!("Model1:{}", model1);
        let model2 = f.generate_with_weight(&dataset, &weight2, k);
        eprintln!("Model2:{}", model2);
        {
            let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
        {
            let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
            eprintln!("{:.4}\t{:.4}", likelihood1, likelihood2);
        }
    }
    use test::Bencher;
    #[bench]
    fn new(b: &mut Bencher) {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let num = 30;
        let template: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..num)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let k = 7;
        b.iter(|| test::black_box(DBGHMM::new(&model1, k)));
    }
    #[bench]
    fn determine(b: &mut Bencher) {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let num = 30;
        let template: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..num)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        let model1 = DBGHMM::new(&model1, k);
        b.iter(|| test::black_box(model1.forward(&template, &DEFAULT_CONFIG)));
    }
    #[bench]
    fn initialize_dbg(b: &mut Bencher) {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let num = 30;
        let template: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..num)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        let model1 = DBGHMM::new(&model1, k);
        b.iter(|| test::black_box(model1.initialize(&template[..k], &DEFAULT_CONFIG)));
    }
    #[bench]
    fn score(b: &mut Bencher) {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let num = 30;
        let template: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..num)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let k = 6;
        let model1 = DBGHMM::new(&model1, k);
        b.iter(|| {
            test::black_box(
                model1
                    .nodes
                    .iter()
                    .map(|e| e.calc_score(&template[..k]))
                    .count(),
            )
        });
    }
    enum Op {
        Match,
        MisMatch,
        Del,
        In,
    }
    impl Op {
        fn weight(&self, p: &Profile) -> f64 {
            match self {
                Op::Match => 1. - p.sub - p.del - p.ins,
                Op::MisMatch => p.sub,
                Op::Del => p.del,
                Op::In => p.ins,
            }
        }
    }
    const OPERATIONS: [Op; 4] = [Op::Match, Op::MisMatch, Op::Del, Op::In];
    fn introduce_randomness<T: rand::Rng>(seq: &[u8], rng: &mut T, p: &Profile) -> Vec<u8> {
        let mut res = vec![];
        let mut remainings: Vec<_> = seq.iter().copied().rev().collect();
        while !remainings.is_empty() {
            match OPERATIONS.choose_weighted(rng, |e| e.weight(p)).unwrap() {
                &Op::Match => res.push(remainings.pop().unwrap()),
                &Op::MisMatch => res.push(choose_base(rng, remainings.pop().unwrap())),
                &Op::In => res.push(random_base(rng)),
                &Op::Del => {
                    remainings.pop().unwrap();
                }
            }
        }
        res
    }
    fn choose_base<T: rand::Rng>(rng: &mut T, base: u8) -> u8 {
        let bases: Vec<u8> = b"ATCG".iter().filter(|&&e| e != base).copied().collect();
        *bases.choose_weighted(rng, |_| 1. / 3.).unwrap()
    }
    fn random_base<T: rand::Rng>(rng: &mut T) -> u8 {
        *b"ATGC".choose_weighted(rng, |_| 1. / 4.).unwrap()
    }
}
