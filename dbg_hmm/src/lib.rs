#![allow(dead_code)]
#![feature(test)]
//! A tiny implementation for de Bruijn graph with Hidden Markov model.
//! Currently, this implementation is minimal. In other words, it exposes only one struct with just two methods:
//! [DeBruijnGraphHiddenMarkovModel] -- Yes, it is too long -- ,[constructor](DeBruijnGraphHiddenMarkovModel::new),
//! and the [forward algorithm](DeBruijnGraphHiddenMarkovModel::forward)
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
const PSEUDO_COUNT: bool = true;
const THR_ON: bool = true;
const THR: f64 = 2.;
const WEIGHT_THR: f64 = 2.0;
const LOW_LIKELIHOOD: f64 = -100000.;
const SCALE: f64 = 3.;
mod find_union;
pub mod gen_sample;
use packed_simd::f64x4 as f64s;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::HashMap;
// This setting is determined by experimentally.
pub const DEFAULT_CONFIG: Config = Config {
    mismatch: 0.03,
    base_freq: [0.25, 0.25, 0.25, 0.25],
    p_match: 0.89,
    p_ins: 0.06,
    p_del: 0.05,
    p_extend_ins: 0.06,
    p_extend_del: 0.05,
    p_del_to_ins: 0.06,
    del_score: -1,
    match_score: 2,
    ins_score: -1,
    mism_score: -1,
};

pub const PACBIO_CONFIG: Config = Config {
    mismatch: 0.0341,
    base_freq: [0.28, 0.22, 0.22, 0.28],
    p_match: 0.9124,
    p_ins: 0.0606,
    p_del: 0.0270,
    p_extend_ins: 0.3014,
    p_extend_del: 0.0695,
    p_del_to_ins: 0.0086,
    del_score: -1,
    match_score: 2,
    ins_score: -1,
    mism_score: -1,
};

pub type DBGHMM = DeBruijnGraphHiddenMarkovModel;
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
            "K:{}\tNodes:{}\tEdges:{}",
            self.k,
            self.nodes.len(),
            edges
        )
    }
}

#[derive(Default)]
pub struct Factory {
    inner: HashMap<Vec<u8>, (f64, usize)>,
    is_safe: Vec<bool>,
    temp_index: Vec<usize>,
    edges: Vec<Vec<usize>>,
}

impl Factory {
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn clear(&mut self) {
        self.inner.clear();
        self.temp_index.clear();
        self.is_safe.clear();
        self.edges.clear();
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
            && self.is_safe.is_empty()
            && self.temp_index.is_empty()
            && self.edges.is_empty()
    }
    pub fn generate(&mut self, dataset: &[Vec<u8>], k: usize) -> DBGHMM {
        // let counter = &mut self.inner;
        dataset.iter().for_each(|seq| {
            seq.windows(k)
                .for_each(|kmer| match self.inner.get_mut(kmer) {
                    Some(x) => x.0 += 1.,
                    None => {
                        self.inner.insert(kmer.to_vec(), (1., std::usize::MAX));
                    }
                })
        });
        let thr = self.calc_thr();
        let mut nodes = Vec::with_capacity(self.len());
        for seq in dataset.iter().filter(|seq| seq.len() > k) {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices(x, &mut nodes, thr, k) {
                    nodes[from].push_edge_with(x[k], to);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let nodes = self.clean_up_nodes(nodes);
        let weight = dataset.len() as f64;
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    pub fn generate_from_ref(&mut self, dataset: &[&[u8]], k: usize) -> DBGHMM {
        assert!(self.is_empty());
        dataset.iter().for_each(|seq| {
            seq.windows(k)
                .for_each(|kmer| match self.inner.get_mut(kmer) {
                    Some(x) => x.0 += 1.,
                    None => {
                        self.inner.insert(kmer.to_vec(), (1., std::usize::MAX));
                    }
                })
        });
        let thr = self.calc_thr();
        let mut nodes = Vec::with_capacity(self.len());
        for seq in dataset {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices(x, &mut nodes, thr, k) {
                    nodes[from].push_edge_with(x[k], to);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let nodes = self.clean_up_nodes(nodes);
        let weight = dataset.len() as f64;
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    pub fn generate_with_weight(&mut self, dataset: &[&[u8]], ws: &[f64], k: usize) -> DBGHMM {
        assert!(self.is_empty());
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(dataset.len() as u64);
        assert_eq!(dataset.len(), ws.len());
        let ep = 0.0001;
        dataset
            .iter()
            .zip(ws.iter())
            .filter(|&(_, &w)| w > ep)
            .for_each(|(seq, w)| {
                seq.windows(k)
                    .for_each(|kmer| match self.inner.get_mut(kmer) {
                        Some(x) => x.0 += w,
                        None => {
                            self.inner.insert(kmer.to_vec(), (*w, std::usize::MAX));
                        }
                    })
            });
        let thr = self.calc_thr_weight();
        let mut nodes = Vec::with_capacity(1_000);
        for (seq, &w) in dataset.iter().zip(ws.iter()).filter(|&(_, &w)| w > ep) {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices_exp(x, &mut nodes, thr, k, &mut rng) {
                    nodes[from].push_edge_with_weight(x[k], to, w);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let weight = ws.iter().sum::<f64>();
        if weight < 1.0001 {
            self.clear();
            let nodes = vec![];
            return DBGHMM { nodes, k, weight };
        }
        let thr = nodes.iter().map(|e| e.tot).sum::<f64>() / nodes.len() as f64 - THR;
        let nodes = match self.clean_up_nodes_exp(nodes, thr, &mut rng) {
            Some(res) => res,
            None => panic!(
                "thr:{},weight:{},raw_kmers:{},wsdump:{:#?}",
                thr,
                weight,
                self.inner.len(),
                ws
            ),
        };
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    fn clean_up_nodes(&mut self, nodes: Vec<Kmer>) -> Vec<Kmer> {
        let mut nodes = self.renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        nodes
    }
    fn clean_up_nodes_exp<R: Rng>(
        &mut self,
        nodes: Vec<Kmer>,
        thr: f64,
        rng: &mut R,
    ) -> Option<Vec<Kmer>> {
        let mut buffer: Vec<_> = Vec::with_capacity(nodes.len());
        let nodes = self.cut_lightweight_loop(nodes, thr, rng, &mut buffer);
        assert!(buffer.is_empty());
        // let nodes = if thr > 0.999 {
        //     let nodes = self.cut_tip(nodes, 2, &mut buffer);
        //     assert!(buffer.is_empty());
        //     self.pick_largest_components(nodes, &mut buffer)
        // } else {
        //     self.pick_largest_components(nodes, &mut buffer)
        // };
        // eprintln!("THR:{}", thr);
        let nodes = self.cut_tip(nodes, 2, &mut buffer);
        assert!(buffer.is_empty());
        let nodes = self.pick_largest_components(nodes, &mut buffer);
        let nodes = match nodes {
            Some(res) => res,
            None => {
                error!("THR:{}", thr);
                return None;
            }
        };
        let mut nodes = self.renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        Some(nodes)
    }
    fn pick_largest_components(
        &mut self,
        mut nodes: Vec<Kmer>,
        buffer: &mut Vec<Kmer>,
    ) -> Option<Vec<Kmer>> {
        assert!(self.temp_index.is_empty());
        assert!(self.is_safe.is_empty());
        let mut fu = find_union::FindUnion::new(nodes.len());
        for (from, n) in nodes.iter().enumerate() {
            for &to in n.edges.iter().filter_map(|e| e.as_ref()) {
                fu.unite(from, to);
            }
        }
        let max_group_and_size = (0..nodes.len())
            .filter_map(|e| {
                let parent = fu.find(e).unwrap();
                if parent == e {
                    fu.size(e).map(|s| (e, s))
                } else {
                    None
                }
            })
            .max_by_key(|&(_, size)| size);
        let (max_group, _) = match max_group_and_size {
            Some(res) => res,
            None => {
                error!("No components:{:?}", nodes);
                return None;
            }
        };
        let mut index = 0;
        for i in 0..nodes.len() {
            self.temp_index.push(index);
            index += (fu.find(i).unwrap() == max_group) as usize;
        }
        let mut index = nodes.len();
        while let Some(mut node) = nodes.pop() {
            index -= 1;
            if fu.find(index).unwrap() == max_group {
                node.remove_if_not(&mut fu, max_group);
                node.rename_by(&self.temp_index);
                buffer.push(node);
            }
        }
        buffer.reverse();
        std::mem::swap(buffer, &mut nodes);
        self.temp_index.clear();
        Some(nodes)
    }
    fn select_supported_node(edges: &[Vec<usize>], nodes: usize, is_safe: &[bool]) -> Vec<bool> {
        let mut is_supported = vec![false; nodes];
        for (from, edges) in edges.iter().enumerate() {
            for &to in edges.iter() {
                is_supported[to] |= is_safe[from];
                is_supported[from] |= is_safe[to];
            }
        }
        is_supported
    }
    fn cut_lightweight_loop<R: Rng>(
        &mut self,
        mut nodes: Vec<Kmer>,
        thr: f64,
        rng: &mut R,
        buffer: &mut Vec<Kmer>,
    ) -> Vec<Kmer> {
        self.is_safe.clear();
        self.temp_index.clear();
        for node in &nodes {
            self.temp_index
                .push(!rng.gen_bool(Self::prob(node.kmer_weight - thr)) as usize)
        }
        (0..nodes.len()).for_each(|_| self.is_safe.push(false));
        for (from, ref node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.is_safe[to] |= self.temp_index[from] == 1;
                self.is_safe[from] |= self.temp_index[to] == 1;
            }
        }
        self.temp_index.clear();
        let mut index = 0;
        for &b in &self.is_safe {
            self.temp_index.push(index);
            index += b as usize;
        }
        assert!(buffer.is_empty());
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            if self.is_safe[idx] {
                node.remove_if_not_supported(&self.is_safe);
                node.rename_by(&self.temp_index);
                buffer.push(node);
            }
        }
        buffer.reverse();
        self.is_safe.clear();
        self.temp_index.clear();
        std::mem::swap(&mut nodes, buffer);
        nodes
    }
    // Cut tips. We can assume that the graph is a
    // connected graph or a tree.
    fn cut_tip(&mut self, mut nodes: Vec<Kmer>, t: usize, buffer: &mut Vec<Kmer>) -> Vec<Kmer> {
        assert!(self.is_safe.is_empty());
        assert!(self.edges.is_empty());
        let pseudo_node = nodes.len();
        (0..nodes.len() + 1).for_each(|_| self.edges.push(vec![]));
        nodes.iter().for_each(|e| self.is_safe.push(e.is_tail));
        for (idx, node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.edges[idx].push(to);
                self.is_safe[idx] |= nodes[to].is_tail;
            }
            if node.is_tail {
                self.edges[idx].push(pseudo_node);
            }
        }
        let is_supported_forward = self.cut_tip_inner(nodes.len() + 1, t);
        self.edges.iter_mut().for_each(|eds| eds.clear());
        nodes
            .iter()
            .zip(self.is_safe.iter_mut())
            .for_each(|(n, b)| *b |= n.is_head);
        for (idx, node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.edges[to].push(idx);
                self.is_safe[to] |= node.is_head;
            }
            if node.is_head {
                self.edges[idx].push(pseudo_node);
            }
        }
        let is_supported_backward = self.cut_tip_inner(nodes.len() + 1, t);
        let ave = nodes.iter().map(|e| e.kmer_weight).sum::<f64>() / nodes.len() as f64;
        nodes
            .iter()
            .zip(self.is_safe.iter_mut())
            .for_each(|(n, b)| *b &= n.kmer_weight > ave / 2.0);
        is_supported_forward
            .into_iter()
            .zip(is_supported_backward)
            .zip(self.is_safe.iter_mut())
            .for_each(|((forward, backward), is_safe)| *is_safe = (forward & backward) | *is_safe);
        assert!(self.temp_index.is_empty());
        let mut idx = 0;
        for &b in &self.is_safe {
            self.temp_index.push(idx);
            idx += b as usize;
        }
        assert!(buffer.is_empty());
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            if self.is_safe[idx] {
                node.remove_if_not_supported(&self.is_safe);
                node.rename_by(&self.temp_index);
                buffer.push(node)
            }
        }
        buffer.reverse();
        std::mem::swap(buffer, &mut nodes);
        self.is_safe.clear();
        self.temp_index.clear();
        self.edges.clear();
        buffer.clear();
        nodes
    }
    // Cut tips. We can assume that the graph is a
    // connected graph or a tree.
    fn cut_tip_inner(&mut self, nodes: usize, t: usize) -> Vec<bool> {
        let edges = &self.edges;
        let t = t as i32;
        let dist_to_root = Self::dist_to_root(edges, nodes);
        let bad_list: Vec<_> = (0..nodes)
            .filter(|&e| edges[e].len() > 1)
            .flat_map(|e| edges[e].iter().filter(|&&to| dist_to_root[to] <= t))
            .copied()
            .collect();
        // If arrived, true.
        let mut is_arrived: Vec<bool> = vec![false; nodes];
        for i in bad_list {
            if is_arrived[i] {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                is_arrived[node] |= true;
                for &to in &edges[node] {
                    if !is_arrived[to] {
                        stack.push(to);
                        continue 'dfs;
                    }
                }
                stack.pop().unwrap();
            }
        }
        is_arrived.iter_mut().for_each(|e| *e = !*e);
        is_arrived
    }
    fn dist_to_root(edges: &[Vec<usize>], nodes: usize) -> Vec<i32> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let mut status = vec![0; nodes];
        let mut dist_to_root: Vec<i32> = vec![-1; nodes];
        for i in 0..nodes {
            if status[i] != 0 {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                if status[node] == 0 {
                    status[node] = 1;
                }
                for &to in &edges[node] {
                    if status[to] == 0 {
                        stack.push(to);
                        continue 'dfs;
                    } else if status[to] == 1 {
                        // Loop detected.
                        dist_to_root[node] = 100000;
                    }
                }
                let last = stack.pop().unwrap();
                dist_to_root[last] = edges[node]
                    .iter()
                    .map(|&to| dist_to_root[to])
                    .max()
                    .unwrap_or(-1)
                    .max(dist_to_root[last])
                    + 1;
                // Deactivate
                status[last] = 2;
            }
        }
        dist_to_root
    }
    fn calc_thr(&self) -> f64 {
        let ave = self.inner.values().map(|e| e.0).sum::<f64>() / self.inner.len() as f64;
        if THR_ON {
            ave - THR
        } else {
            0.
        }
    }
    fn calc_thr_weight(&self) -> f64 {
        let (sum, denom) =
            self.inner.values().fold(
                (0., 0.),
                |(x, y), &(w, _)| if w >= 1. { (x + w, y + 1.) } else { (x, y) },
            );
        let ave = sum / denom;
        if THR_ON {
            // y = x - 1 - epsilon
            ave - THR
        } else {
            0.
        }
    }
    fn get_idx(&mut self, kmer: &[u8], nodes: &mut Vec<Kmer>, thr: f64) -> Option<usize> {
        match self.inner.get_mut(kmer) {
            Some(x) if x.0 <= thr => None,
            Some(x) if x.1 != std::usize::MAX => Some(x.1),
            Some(x) => {
                x.1 = nodes.len();
                nodes.push(Kmer::new(kmer, x.0));
                Some(nodes.len() - 1)
            }
            _ => unreachable!(),
        }
    }
    fn prob(x: f64) -> f64 {
        // if x.is_sign_positive(){
        //     0.
        // }else{
        ((SCALE * x).exp() + 1.).recip()
        // }
    }
    fn get_idx_exp<R: Rng>(
        &mut self,
        kmer: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        rng: &mut R,
    ) -> Option<usize> {
        match self.inner.get_mut(kmer) {
            Some(x) if rng.gen_bool(Self::prob(x.0 - thr)) => None,
            Some(x) if x.1 != std::usize::MAX => Some(x.1),
            Some(x) => {
                x.1 = nodes.len();
                nodes.push(Kmer::new(kmer, x.0));
                Some(nodes.len() - 1)
            }
            _ => unreachable!(),
        }
    }

    fn get_indices(
        &mut self,
        x: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        k: usize,
    ) -> Option<(usize, usize)> {
        let (prev, after) = (&x[..k], &x[1..]);
        let from = self.get_idx(prev, nodes, thr)?;
        let to = self.get_idx(after, nodes, thr)?;
        Some((from, to))
    }
    fn get_indices_exp<R: Rng>(
        &mut self,
        x: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        k: usize,
        rng: &mut R,
    ) -> Option<(usize, usize)> {
        let (prev, next) = (&x[..k], &x[1..]);
        let from = self.get_idx_exp(prev, nodes, thr, rng)?;
        let next = self.get_idx_exp(next, nodes, thr, rng)?;
        Some((from, next))
    }

    // Return topological sorted order. If there's a cycle,
    // it neglect the back-edge, and proceed with raise error messages to stderr(DEBUG MODE).
    fn topological_sort(nodes: &[Kmer]) -> Vec<usize> {
        let mut edges: Vec<Vec<_>> = vec![vec![]; nodes.len()];
        for (idx, node) in nodes.iter().enumerate() {
            for to in node.edges.iter().filter_map(|e| e.as_ref()) {
                edges[idx].push(*to);
            }
        }
        match Self::topological_sort_inner(&edges, nodes.len()) {
            Ok(res) => res,
            Err(res) => res,
        }
    }
    // Topological sorting.
    fn topological_sort_inner(
        edges: &[Vec<usize>],
        nodes: usize,
    ) -> Result<Vec<usize>, Vec<usize>> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let mut status = vec![0; nodes];
        let mut order = vec![];
        let mut is_dag = true;
        for i in 0..nodes {
            if status[i] != 0 {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                if status[node] == 0 {
                    // preorder
                    status[node] = 1;
                }
                for &to in &edges[node] {
                    if status[to] == 0 {
                        stack.push(to);
                        continue 'dfs;
                    } else if status[to] == 1 {
                        is_dag = false;
                    }
                }
                // No-op
                let last = stack.pop().unwrap();
                order.push(last);
                // Deactivate
                status[last] = 2;
            }
        }
        order.reverse();
        if is_dag {
            Ok(order)
        } else {
            Err(order)
        }
    }
    fn renaming_nodes(&mut self, mut nodes: Vec<Kmer>) -> Vec<Kmer> {
        assert!(self.temp_index.is_empty());
        (0..nodes.len()).for_each(|_| self.temp_index.push(0));
        let order = Self::topological_sort(&nodes);
        for (order, v) in order.into_iter().enumerate() {
            self.temp_index[v] = order;
        }
        let mut result = vec![None; nodes.len()];
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            node.rename_by(&self.temp_index);
            result[self.temp_index[idx]] = Some(node);
        }
        assert!(result.iter().all(|e| e.is_some()));
        assert!(nodes.is_empty());
        self.temp_index.clear();
        nodes.extend(result.into_iter().filter_map(|e| e));
        nodes
    }
    pub fn new() -> Self {
        let inner = std::collections::HashMap::new();
        let is_safe = vec![];
        let temp_index = vec![];
        let edges = vec![];
        Self {
            inner,
            is_safe,
            temp_index,
            edges,
        }
    }
}

impl DeBruijnGraphHiddenMarkovModel {
    pub fn new(dataset: &[Vec<u8>], k: usize) -> Self {
        let mut f = Factory::new();
        f.generate(dataset, k)
    }
    pub fn new_from_ref(dataset: &[&[u8]], k: usize) -> Self {
        let mut f = Factory::new();
        f.generate_from_ref(dataset, k)
    }
    // Calc hat. Return (c,d)
    #[cfg(not(target_feature = "sse"))]
    fn update(&self, updates: &mut [Node], prev: &[Node], base: u8, config: &Config) -> (f64, f64) {
        for (idx, (from, node)) in self
            .nodes
            .iter()
            .zip(prev.iter())
            .enumerate()
            .filter(|(_, (_, node))| !node.is_zero())
        {
            updates[idx].add_ins(node.insertion(&config) * from.insertion(base));
            let prob = node.match_succeed(config);
            from
                .edges
                .iter()
                .enumerate()
                .for_each(|(i,e)|
                          // (from->to,base i edge)
                          if let Some(to) = e {
                              updates[*to].add_mat(prob * from.to(i) * self.nodes[*to].prob(base, config))
                          });
        }
        let d = Node::sum(&updates).recip();
        // Updates -> tilde
        updates.iter_mut().for_each(|e| *e = *e * d);
        // Update `Del` states.
        for (idx, node) in self.nodes.iter().enumerate() {
            if updates[idx].is_nonzero() {
                let pushing_weight = updates[idx].push(&config);
                for to in node.edges.iter() {
                    match to {
                        Some(to) => updates[*to].add_del(pushing_weight),
                        None => {}
                    }
                }
            }
        }
        let c = Node::sum(&updates).recip();
        updates.iter_mut().for_each(|e| *e = *e * c);
        (c, d)
    }

    #[cfg(target_feature = "sse")]
    fn is_non_zero(idx: usize, xs: &[f64]) -> bool {
        xs[idx * 4..(idx + 1) * 4].iter().any(|&e| e > 0.00001)
    }
    #[cfg(target_feature = "sse")]
    fn sum(xs: &[f64]) -> f64 {
        assert!(xs.len() % 4 == 0);
        xs.chunks_exact(f64s::lanes())
            .map(f64s::from_slice_unaligned)
            .sum::<f64s>()
            .sum()
    }
    #[cfg(target_feature = "sse")]
    fn mul(xs: &mut [f64], y: f64) {
        assert!(xs.len() % 4 == 0);
        let ys = f64s::splat(y);
        xs.chunks_exact_mut(f64s::lanes()).for_each(|xs| {
            let packed = f64s::from_slice_unaligned(xs) * ys;
            packed.write_to_slice_unaligned(xs);
        });
    }
    #[cfg(target_feature = "sse")]
    fn update(&self, updates: &mut [f64], prev: &[f64], base: u8, config: &Config) -> (f64, f64) {
        // Alignemnt:[mat,ins,del,0., mat,ins,del,0., mat,....,del, 0.,]
        for (idx, from) in self
            .nodes
            .iter()
            .enumerate()
            .filter(|&(idx, _)| Self::is_non_zero(idx, prev))
        {
            let node = 4 * idx;
            let prob = prev[node] * config.p_match
                + prev[node + 1] * (1. - config.p_extend_ins)
                + prev[node + 2] * (1. - config.p_extend_del - config.p_del_to_ins);
            let ins = prev[node] * config.p_ins
                + prev[node + 1] * config.p_extend_ins
                + prev[node + 2] * config.p_del_to_ins;
            updates[node + 1] += ins * from.insertion(base);
            from
                .edges
                .iter()
                .enumerate()
                .for_each(|(i,e)|
                          // (from->to,base i edge)
                          if let Some(to) = *e {
                              updates[4*to]+= prob * from.to(i) * self.nodes[to].prob(base, config)
                          });
        }
        let d = Self::sum(updates).recip();
        // Updates -> tilde
        Self::mul(updates, d);
        // Update `Del` states.
        for (idx, node) in self.nodes.iter().enumerate() {
            if Self::is_non_zero(idx, updates) {
                let from = &updates[4 * idx..4 * (idx + 1)];
                let pushing_weight: f64 = from[0] * config.p_del + from[2] * config.p_extend_del;
                for to in node.edges.iter().map(|e| e.as_ref()) {
                    if let Some(&to) = to {
                        updates[4 * to + 2] += pushing_weight;
                    }
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
    fn initialize(&self, tip: &[u8], config: &Config) -> (f64, f64, Vec<Node>) {
        // assert!(self.nodes.iter().any(|e| e.is_head));
        let alignments: Vec<_> = self.nodes.iter().map(|e| e.calc_score(tip)).collect();
        let initial_prob: Vec<_> = vec![(self.nodes.len() as f64).recip(); self.nodes.len()];
        let aln_scores: Vec<_> = alignments
            .iter()
            .map(|&(_, mis, del, ins)| {
                config.mismatch.ln() * (mis as f64)
                    + config.p_del.ln() * (del as f64)
                    + config.p_ins.ln() * (ins as f64)
            })
            .zip(initial_prob)
            .map(|(x, y)| x + y.ln())
            .collect();
        let minus_ln_d = Self::logsumexp(&aln_scores);
        let tilde = aln_scores.into_iter().map(|x| (x - minus_ln_d).exp());
        let hat: Vec<_> = tilde.map(|init| Node::new(init, 0., 0.)).collect();
        let minus_ln_c = 0.;
        (minus_ln_c, minus_ln_d, hat)
    }
    /// Return the total weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    // This returns log p(obs|model) = \sum - log c_t.
    #[cfg(not(target_feature = "sse"))]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        if self.weight < WEIGHT_THR {
            return LOW_LIKELIHOOD;
        }
        let (mut cs, mut ds, mut prev) = self.initialize(&obs[..self.k], config);
        let mut updated = vec![Node::new(0., 0., 0.); self.nodes.len()];
        for &base in &obs[self.k..] {
            updated.iter_mut().for_each(Node::clear);
            let (c, d) = self.update(&mut updated, &prev, base, config);
            cs -= c.ln();
            ds -= d.ln();
            std::mem::swap(&mut prev, &mut updated);
        }
        assert!(cs + ds < 0.);
        cs + ds
    }
    #[cfg(target_feature = "sse")]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        if self.weight < WEIGHT_THR {
            return LOW_LIKELIHOOD;
        }
        let (mut cs, mut ds, prev) = self.initialize(&obs[..self.k], config);
        // Alignemnts: [mat, ins, del, 0, mat, ins, del, 0, ....]
        let mut prev: Vec<f64> = prev.iter().flat_map(|e| vec![e.mat, 0., 0., 0.]).collect();
        let mut updated = vec![0.; self.node_num() * 4];
        for &base in &obs[self.k..] {
            updated.iter_mut().for_each(|e| *e = 0.);
            let (c, d) = self.update(&mut updated, &prev, base, config);
            cs -= c.ln();
            ds -= d.ln();
            std::mem::swap(&mut prev, &mut updated);
        }
        assert!(cs + ds < 0.);
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
}
#[derive(Clone, Copy)]
struct Node {
    mat: f64,
    del: f64,
    ins: f64,
}

impl std::fmt::Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:.4}\t{:.4}\t{:.4}", self.mat, self.del, self.ins)
    }
}

impl std::ops::Mul<f64> for Node {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let Node { mat, del, ins } = self;
        Node {
            mat: mat * rhs,
            del: del * rhs,
            ins: ins * rhs,
        }
    }
}

impl Node {
    fn sum(xs: &[Node]) -> f64 {
        xs.iter().map(Node::fold).sum::<f64>()
    }
    fn new(mat: f64, del: f64, ins: f64) -> Self {
        Self { mat, del, ins }
    }
    fn is_zero(&self) -> bool {
        self.mat < 0.00001 && self.del < 0.00001 && self.ins < 0.00001
    }
    fn is_nonzero(&self) -> bool {
        self.mat > 0.00001 || self.del > 0.00001 || self.ins > 0.00001
    }
    fn clear(&mut self) {
        self.mat = 0.;
        self.del = 0.;
        self.ins = 0.;
    }
    fn match_succeed(&self, config: &Config) -> f64 {
        self.mat * config.p_match
            + self.del * (1. - config.p_extend_del - config.p_del_to_ins)
            + self.ins * (1. - config.p_extend_ins)
    }
    fn insertion(&self, config: &Config) -> f64 {
        (self.mat * config.p_ins + self.ins * config.p_extend_ins + self.del * config.p_del_to_ins)
    }
    fn push(&self, config: &Config) -> f64 {
        self.mat * config.p_del + self.del * config.p_extend_del
    }
    fn add_ins(&mut self, x: f64) {
        self.ins += x;
    }
    fn add_del(&mut self, x: f64) {
        self.del += x;
    }
    fn add_mat(&mut self, x: f64) {
        self.mat += x;
    }
    fn fold(&self) -> f64 {
        self.mat + self.del + self.ins
    }
}

/// A configure struct.
#[derive(Debug, Clone, Default)]
pub struct Config {
    /// Mismatch probability at given position. # mismatch/(#mism + #match)
    pub mismatch: f64,
    /// The base composition of the read. A,C,G, and T.
    pub base_freq: [f64; 4],
    /// probability of matching at match states.
    pub p_match: f64,
    /// probability of starting insertion at match states.
    pub p_ins: f64,
    /// probability of starting deletion at match states.
    pub p_del: f64,
    /// probability of extending insertion at insertion state.
    pub p_extend_ins: f64,
    /// same for deletion
    pub p_extend_del: f64,
    /// probability of jump deletion state to insertion state.
    pub p_del_to_ins: f64,
    /// match score
    pub match_score: i32,
    /// mismatch score
    pub mism_score: i32,
    /// deletion score
    pub del_score: i32,
    /// insertion score
    pub ins_score: i32,
}

impl Config {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mismatch: f64,
        p_match: f64,
        p_ins: f64,
        p_del: f64,
        p_extend_ins: f64,
        p_extend_del: f64,
        p_del_to_ins: f64,
        match_score: i32,
        ins_score: i32,
        del_score: i32,
        mism_score: i32,
        base_freq: [f64; 4],
    ) -> Self {
        Self {
            mism_score,
            match_score,
            del_score,
            ins_score,
            mismatch,
            base_freq,
            p_match,
            p_ins,
            p_del,
            p_extend_ins,
            p_extend_del,
            p_del_to_ins,
        }
    }
}

#[derive(Clone, Default)]
struct Kmer {
    kmer: Vec<u8>,
    last: u8,
    weight: [f64; 4],
    transition: [f64; 4],
    // Total number of *Outdegree*
    tot: f64,
    // The location to the edges with the label of A,C,G,and T.
    // If there is no edges, None
    edges: [Option<usize>; 4],
    // Weight of this kmer.
    kmer_weight: f64,
    // Whether this is the end of unit.
    is_tail: bool,
    is_head: bool,
}

impl std::fmt::Debug for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Kmer:{}", String::from_utf8_lossy(&self.kmer))?;
        writeln!(f, "KmerWeight:{:.3}", self.kmer_weight)?;
        writeln!(f, "Last:{}", self.last as char)?;
        writeln!(f, "Weight:{:?}", self.weight)?;
        writeln!(f, "Transition:{:?}", self.transition)?;
        writeln!(f, "tot:{}", self.tot)?;
        writeln!(f, "is_tail:{}", self.is_tail)?;
        writeln!(f, "is_head:{}", self.is_head)?;
        for (i, to) in self
            .edges
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| e.map(|to| (idx, to)))
        {
            writeln!(f, "->{}({})", to, i)?;
        }
        Ok(())
    }
}

impl Kmer {
    fn new(x: &[u8], w: f64) -> Self {
        let kmer = x.to_vec();
        let kmer_weight = w;
        let last = *kmer.last().unwrap();
        // Prior
        let weight = if PSEUDO_COUNT { [1.; 4] } else { [0.; 4] };
        let tot = 0.;
        let transition = [0f64; 4];
        let edges = [None; 4];
        let is_tail = false;
        let is_head = false;
        Self {
            kmer,
            kmer_weight,
            last,
            weight,
            tot,
            edges,
            transition,
            is_tail,
            is_head,
        }
    }
    fn finalize(&mut self) {
        let tot_for_weight = if PSEUDO_COUNT {
            self.tot + 4.
        } else {
            self.tot
        };
        if self.tot > 0. {
            for i in 0..4 {
                self.weight[i] /= tot_for_weight;
                self.transition[i] /= self.tot;
            }
        // assert!((1. - self.transition.iter().sum::<f64>()).abs() < 0.0001);
        // assert!((1. - self.weight.iter().sum::<f64>()).abs() < 0.0001);
        } else {
            self.weight = [0.25; 4];
        }
    }
    // renaming all the edges by `map`
    fn rename_by(&mut self, map: &[usize]) {
        for edge in self.edges.iter_mut() {
            if let Some(res) = edge.as_mut() {
                *res = map[*res];
            }
        }
    }
    // Remove all the edges to nodes not in the maximum group.
    fn remove_if_not(&mut self, fu: &mut find_union::FindUnion, mg: usize) {
        for i in 0..4 {
            if let Some(res) = self.edges[i] {
                if fu.find(res).unwrap() != mg {
                    self.edges[i] = None;
                    self.tot -= self.transition[i];
                    self.weight[i] -= self.transition[i];
                    self.transition[i] = 0.;
                }
            }
        }
    }
    // Remove all the edges to unsuppoeted nodes.
    fn remove_if_not_supported(&mut self, is_supported: &[bool]) {
        for i in 0..4 {
            if let Some(res) = self.edges[i] {
                if !is_supported[res] {
                    self.edges[i] = None;
                    self.tot -= self.transition[i];
                    self.weight[i] -= self.transition[i];
                    self.transition[i] = 0.;
                }
            }
        }
    }

    fn push_edge_with(&mut self, base: u8, to: usize) {
        let i = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        };
        match self.edges.get_mut(i) {
            Some(Some(target)) => assert_eq!(*target, to),
            Some(x) => *x = Some(to),
            _ => unreachable!(),
        };
        self.tot += 1.;
        self.weight[i] += 1.;
        self.transition[i] += 1.;
    }
    fn push_edge_with_weight(&mut self, base: u8, to: usize, w: f64) {
        let i = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        };
        match self.edges.get_mut(i) {
            Some(Some(target)) => assert_eq!(*target, to),
            Some(x) => *x = Some(to),
            _ => unreachable!(),
        };
        self.tot += w;
        self.weight[i] += w;
        self.transition[i] += w;
    }
    // return (Score(self.kmer,tip), match, mism,del, ins)
    fn calc_score(&self, tip: &[u8]) -> (u8, u8, u8, u8) {
        let aln = edlib_sys::global(&self.kmer, tip);
        aln.into_iter()
            .fold((0, 0, 0, 0), |(mat, mis, del, ins), x| match x {
                0 => (mat + 1, mis, del, ins),
                1 => (mat, mis, del, ins + 1),
                2 => (mat, mis, del + 1, ins),
                3 => (mat, mis + 1, del, ins),
                _ => unreachable!(),
            })
    }
    // return P(idx|self)
    #[inline]
    fn to(&self, idx: usize) -> f64 {
        self.transition[idx]
    }
    #[inline]
    fn prob(&self, base: u8, config: &Config) -> f64 {
        if self.last == base {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        }
    }
    // return P_I(base|self)
    fn insertion(&self, base: u8) -> f64 {
        match base {
            b'A' => self.weight[0],
            b'C' => self.weight[1],
            b'G' => self.weight[2],
            b'T' => self.weight[3],
            _ => 0.,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::Xoroshiro128StarStar;
    #[test]
    fn works() {}
    #[test]
    fn tip_cut() {
        let edges = vec![
            vec![1, 2],
            vec![],
            vec![],
            vec![4],
            vec![0, 5],
            vec![9],
            vec![7],
            vec![8],
            vec![11],
            vec![6, 7, 10],
            vec![],
            vec![],
        ];
        let nodes = 12;
        assert_eq!(edges.len(), 12);
        let dist_to_root = Factory::dist_to_root(&edges, nodes);
        let answer = vec![1, 0, 0, 7, 6, 5, 3, 2, 1, 4, 0, 0];
        assert_eq!(dist_to_root, answer);
        // let is_supported = Factory::cut_tip_inner(&edges, nodes, 1 );
        // let answer = vec![
        //     false, false, false, true, true, true, true, true, true, true, false, true,
        // ];
        // assert_eq!(is_supported, answer);
        let edges = vec![
            vec![1],
            vec![2, 3],
            vec![],
            vec![4],
            vec![],
            vec![6],
            vec![7],
            vec![0, 8],
            vec![5],
        ];
        let nodes = 9;
        assert_eq!(edges.len(), nodes);
        let dist_to_root = Factory::dist_to_root(&edges, nodes);
        let answer = vec![3, 2, 0, 1, 0];
        eprintln!("{:?}", dist_to_root);
        for i in 0..5 {
            assert_eq!(dist_to_root[i], answer[i]);
        }
        for i in 5..nodes {
            assert!(dist_to_root[i] > 100);
        }
        // let is_supported = Factory::cut_tip_inner(&edges, nodes, 1);
        // let answer = vec![true, true, false, false, false, true, true, true, true];
        // assert_eq!(answer, is_supported);
    }
    #[test]
    fn select_supported_test() {
        let edges = vec![
            vec![1],
            vec![2],
            vec![3, 7],
            vec![],
            vec![5],
            vec![],
            vec![7],
            vec![4],
            vec![6],
            vec![8, 0],
            vec![9],
        ];
        let nodes = 11;
        let is_safe = vec![
            false, true, true, false, true, true, false, false, false, true, true,
        ];
        let supported = Factory::select_supported_node(&edges, nodes, &is_safe);
        let mut answer = vec![true; nodes];
        answer[6] = false;
        assert_eq!(supported, answer);
    }
    #[test]
    fn topological_sorting() {
        let edges = vec![
            vec![1],
            vec![2],
            vec![3],
            vec![4],
            vec![],
            vec![6],
            vec![2],
            vec![8],
            vec![9],
            vec![3],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap();
        assert_eq!(order, vec![7, 8, 9, 5, 6, 0, 1, 2, 3, 4]);
        let edges = vec![
            vec![1],
            vec![2],
            vec![3],
            vec![4, 8],
            vec![],
            vec![6],
            vec![2],
            vec![8],
            vec![9],
            vec![3],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap_err();
        assert_eq!(order, vec![7, 5, 6, 0, 1, 2, 3, 8, 9, 4]);
        let edges = vec![
            vec![1],
            vec![2, 7],
            vec![3],
            vec![4],
            vec![5],
            vec![6],
            vec![],
            vec![8],
            vec![5],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap();
        assert_eq!(order, vec![0, 1, 7, 8, 2, 3, 4, 5, 6]);
    }
    #[test]
    fn kmer() {
        let mut kmer = Kmer::new(b"ACTCGTA", 0.);
        assert_eq!(kmer.last, b'A');
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'C', 34);
        kmer.push_edge_with(b'G', 32);
        kmer.push_edge_with(b'A', 12);
        if PSEUDO_COUNT {
            assert_eq!(kmer.weight, [1. + 1. + 1., 1. + 1., 1. + 1., 1.]);
        } else {
            assert_eq!(kmer.weight, [1. + 1., 1., 1., 0.]);
        }
        assert_eq!(kmer.edges, [Some(12), Some(34), Some(32), None]);
    }
    #[test]
    #[should_panic]
    fn kmer2() {
        let mut kmer = Kmer::new(b"ACTCGTA", 1.);
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'A', 34);
    }
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
        sub: 0.02,
        del: 0.065,
        ins: 0.065,
    };
    const PAC_BIO: Profile = Profile {
        sub: 0.03,
        del: 0.04,
        ins: 0.07,
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
            let model1: Vec<Vec<_>> = (0..num)
                .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
                .collect();
            let weight = vec![1.; num];
            let kmers: std::collections::HashSet<Vec<u8>> = model1
                .iter()
                .flat_map(|e| e.windows(k))
                .map(|e| e.to_vec())
                .collect();
            let model1: Vec<_> = model1.iter().map(|e| e.as_slice()).collect();
            let mut f = Factory::new();
            let model1 = f.generate_with_weight(&model1, &weight, k);
            eprintln!("{}", model1);
            eprintln!("{}", String::from_utf8_lossy(&template1));
            eprintln!("{}", num);
            let num = template1
                .windows(k)
                .filter(|&kmer| !model1.nodes.iter().any(|node| node.kmer == kmer))
                .inspect(|kmer| eprintln!("{}", String::from_utf8_lossy(kmer)))
                .count();
            eprintln!("{}", num);
            for kmer in template1.windows(k).filter(|&kmer| kmers.contains(kmer)) {
                assert!(
                    model1.nodes.iter().any(|node| node.kmer == kmer) // "{:?}",
                                                                      // model1
                );
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
    fn single_error_test() {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1234565);
        let coverage = 50;
        let results: Vec<_> = (7..coverage)
            .step_by(2)
            .map(|cov| {
                let template1: Vec<_> = (0..150)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 1, 0, 0);
                let sub = check(&template1, &template2, &mut rng, cov);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 1, 0);
                let del = check(&template1, &template2, &mut rng, cov);
                let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 0, 1);
                let ins = check(&template1, &template2, &mut rng, cov);
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
        assert!(false);
    }
    fn check<R: rand::Rng>(t1: &[u8], t2: &[u8], rng: &mut R, cov: usize) -> usize {
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
        let mut f = Factory::new();
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
