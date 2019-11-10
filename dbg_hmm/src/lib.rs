#![feature(test)]
//! A tiny implementation for de Bruijn graph with Hidden Markov model.
//! Currently, this implementation is minimal. In other words, it exposes only one struct with just two methods:
//! [DeBruijnGraphHiddenMarkovModel] -- Yes, it is too long -- ,[constructor](DeBruijnGraphHiddenMarkovModel::new),
//! and the [forward algorithm](DeBruijnGraphHiddenMarkovModel::forward)
//! to calculate the probability this graph would generate the given observation.
//! As a shorthand for the vary long name, I also supply [DBGHMM] as a alias for [DeBruijnGraphHiddenMarkovModel].
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate test;
#[macro_use]
extern crate log;
extern crate env_logger;
// Whether or not to use 'pseudo count' in the out-dgree.
const PSEUDO_COUNT: bool = true;
pub mod gen_sample;
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
    p_del_to_ins: 0.86,
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
}

pub struct Factory {
    inner: std::collections::HashMap<Vec<u8>, usize>,
}
impl Factory {
    pub fn generate(&mut self, dataset: &[Vec<u8>], k: usize) -> DBGHMM {
        let indexer = &mut self.inner;
        let mut nodes = vec![];
        for seq in dataset {
            for x in seq.windows(k + 1) {
                let from = Self::push(indexer, &mut nodes, &x[..k]);
                let to = Self::push(indexer, &mut nodes, &x[1..]);
                nodes[from].push_edge_with(x[k], to);
            }
        }
        // Node Index -> Topological order.
        let mut nodes = Self::renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        self.inner.clear();
        DBGHMM { nodes, k }
    }
    pub fn generate_from_ref(&mut self, dataset: &[&[u8]], k: usize) -> DBGHMM {
        let indexer = &mut self.inner;
        let mut nodes = vec![];
        for seq in dataset {
            for x in seq.windows(k + 1) {
                let from = Self::push(indexer, &mut nodes, &x[..k]);
                let to = Self::push(indexer, &mut nodes, &x[1..]);
                nodes[from].push_edge_with(x[k], to);
            }
        }
        let mut nodes = Self::renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        self.inner.clear();
        DBGHMM { nodes, k }
    }
    fn push(hm: &mut HashMap<Vec<u8>, usize>, nodes: &mut Vec<Kmer>, kmer: &[u8]) -> usize {
        if !hm.contains_key(kmer) {
            hm.insert(kmer.to_vec(), nodes.len());
            nodes.push(Kmer::new(kmer));
            nodes.len() - 1
        } else {
            hm[kmer]
        }
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
            Err(res) => {
                debug!("The graph is cyclic.");
                debug!("{:?}", nodes);
                res
            }
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
    fn renaming_nodes(mut nodes: Vec<Kmer>) -> Vec<Kmer> {
        let order = Self::topological_sort(&nodes);
        let rev_index = {
            let mut ri = vec![0; nodes.len()];
            for (order, v) in order.into_iter().enumerate() {
                ri[v] = order;
            }
            ri
        };
        let mut result = vec![None; nodes.len()];
        nodes.iter_mut().for_each(|e| e.rename_by(&rev_index));
        for (idx, n) in nodes.into_iter().enumerate() {
            result[rev_index[idx]] = Some(n);
        }
        assert!(result.iter().all(|e| e.is_some()));
        result.into_iter().filter_map(|e| e).collect()
    }
    pub fn new() -> Self {
        let inner = std::collections::HashMap::new();
        Self { inner }
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
    fn update(&self, updates: &mut [Node], prev: &[Node], base: u8, config: &Config) -> (f64, f64) {
        // let st = std::time::Instant::now();
        for (idx, (from, node)) in self.nodes.iter().zip(prev.iter()).enumerate() {
            if node.is_zero() {
                continue;
            }
            let prob = node.match_succeed(config);
            // Update `Mat` states. updates -> double dots
            for (i, to) in from
                .edges
                .iter()
                .enumerate()
                .filter_map(|(i, e)| e.map(|n| (i, n)))
            {
                // (from->to,base i edge)
                updates[to].add_mat(prob * from.to(i) * self.nodes[to].prob(base, config));
            }
            // Update `Ins` states. updates -> double dots
            updates[idx].add_ins(node.insertion(&config) * self.nodes[idx].insertion(base));
        }

        // Calculate D.
        let d = 1. / updates.iter().map(|e| e.mat + e.ins).sum::<f64>();
        // Updates -> tilde
        updates.iter_mut().for_each(|e| *e = *e * d);
        // Update `Del` states.
        // Up to thi line, almost free.
        assert!(updates.iter().all(|e| e.del == 0.));
        let non_zero_nodes: Vec<_> = self
            .nodes
            .iter()
            .zip(updates.iter())
            .filter(|(_, from)| !from.is_zero())
            .map(|(n, from)| (n, from.push(&config)))
            .collect();
        for (node, pushing_weight) in non_zero_nodes {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                updates[to].add_del(pushing_weight);
            }
        }
        // calculate c.
        let c = 1. / updates.iter().map(Node::fold).sum::<f64>();
        updates.iter_mut().for_each(|e| *e = *e * c);
        // let st2 = std::time::Instant::now();
        // eprintln!("{:?}\t{:?}", st2 - st , st3 - st2);
        (c, d)
    }
    fn initialize(&self, tip: &[u8], config: &Config) -> (f64, f64, Vec<Node>) {
        //debug!("Query:{}", String::from_utf8_lossy(tip));
        let initial_prob: Vec<_> = self
            .nodes
            .iter()
            .map(|e| (e.calc_score(tip, config) as f64).exp())
            .collect();
        let s = initial_prob.iter().sum::<f64>();
        let initial_prob = initial_prob.into_iter().map(|e| e / s);
        let last = tip[tip.len() - 1];
        let double_dots: Vec<_> = self
            .nodes
            .iter()
            .zip(initial_prob)
            .map(|(node, init)| node.prob(last, config) * init)
            .collect();
        let d = 1. / double_dots.iter().map(|e| e).sum::<f64>();
        let tilde = double_dots.into_iter().map(|e| e * d);
        // let c = 1. / tilde.iter().sum::<f64>(); Actually, the below is true.
        let c = 1.;
        let hat: Vec<_> = tilde.map(|init| Node::new(init, 0., 0.)).collect();
        (c, d, hat)
    }
    // This returns log p(obs|model) = \sum - log c_t.
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        //debug!("Nodes:{:?}", self.nodes);
        let mut cs = Vec::with_capacity(obs.len() + 1);
        let mut ds = Vec::with_capacity(obs.len() + 1);
        // let's initialize the first array.
        let (c, d, mut prev) = self.initialize(&obs[..self.k], config);
        cs.push(c);
        ds.push(d);
        let mut updated = vec![Node::new(0., 0., 0.); self.nodes.len()];
        // let mut idx = self.k;
        for &base in &obs[self.k..] {
            // debug!("{} nodes, weight sum:{}", self.nodes.len(), csum);
            // idx += 1;
            // debug!("{}-th base {}", idx, base as char);
            updated.iter_mut().for_each(Node::clear);
            // Calculate hat.
            // Calculate c. Calculate d.
            let (c, d) = self.update(&mut updated, &prev, base, config);
            cs.push(c);
            ds.push(d);
            // Flip previous and updated.
            std::mem::swap(&mut prev, &mut updated);
        }
        -cs.into_iter().map(|e| e.ln()).sum::<f64>() - ds.into_iter().map(|e| e.ln()).sum::<f64>()
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
    fn new(mat: f64, del: f64, ins: f64) -> Self {
        Self { mat, del, ins }
    }
    fn is_zero(&self) -> bool {
        self.mat < 0.00001 && self.del < 0.00001 && self.ins < 0.00001
    }
    #[inline]
    fn clear(&mut self) {
        self.mat = 0.;
        self.del = 0.;
        self.ins = 0.;
    }
    #[inline]
    fn match_succeed(&self, config: &Config) -> f64 {
        self.mat * config.p_match
            + self.del * (1. - config.p_extend_del - config.p_del_to_ins)
            + self.ins * (1. - config.p_extend_ins)
    }
    fn insertion(&self, config: &Config) -> f64 {
        (self.mat * config.p_ins + self.ins * config.p_extend_ins + self.del * config.p_del_to_ins)
    }
    fn push(&self, config: &Config) -> f64 {
        self.del * config.p_extend_del + self.mat * config.p_del
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
    mismatch: f64,
    /// The base composition of the read. A,C,G, and T.
    base_freq: [f64; 4],
    /// probability of matching at match states.
    p_match: f64,
    /// probability of starting insertion at match states.
    p_ins: f64,
    /// probability of starting deletion at match states.
    p_del: f64,
    /// probability of extending insertion at insertion state.
    p_extend_ins: f64,
    /// same for deletion
    p_extend_del: f64,
    /// probability of jump deletion state to insertion state.
    p_del_to_ins: f64,
    /// match score
    match_score: i32,
    /// mismatch score
    mism_score: i32,
    /// deletion score
    del_score: i32,
    /// insertion score
    ins_score: i32,
}

impl Config {
    pub fn new(
        mismatch: f64,
        base_freq: [f64; 4],
        p_match: f64,
        p_ins: f64,
        p_del: f64,
        p_extend_ins: f64,
        p_extend_del: f64,
        p_del_to_ins: f64,
        del_score: i32,
        match_score: i32,
        ins_score: i32,
        mism_score: i32,
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
    weight: [u16; 4],
    tot: u16,
    probs: [f64; 4],
    // The location to the edges with the label of A,C,G,and T.
    // If there is no edges, None
    edges: [Option<usize>; 4],
}

impl std::fmt::Debug for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Kmer:{}", String::from_utf8_lossy(&self.kmer))?;
        writeln!(f, "Last:{}", self.last as char)?;
        writeln!(f, "Weight:{:?}", self.weight)?;
        writeln!(f, "tot:{}", self.tot)?;
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
    fn new(x: &[u8]) -> Self {
        let kmer = x.to_vec();
        let last = *kmer.last().unwrap();
        // Prior
        let (weight, tot) = if PSEUDO_COUNT {
            ([1; 4], 4)
        } else {
            ([0; 4], 0)
        };
        let edges = [None; 4];
        let probs = [0.25; 4];
        Self {
            kmer,
            last,
            weight,
            tot,
            edges,
            probs,
        }
    }
    #[allow(dead_code)]
    fn polish(self) -> Self {
        let mut edges = [None; 4];
        let mut weight = self.weight;
        let mut tot = self.tot;
        for i in 0..4 {
            edges[i] = match self.edges[i] {
                Some(res) if weight[i] > 1 => Some(res),
                Some(_) => {
                    tot -= weight[i];
                    weight[i] = 0;
                    None
                }
                None => None,
            };
        }
        let mut probs = self.probs;
        for i in 0..4 {
            probs[i] = weight[i] as f64 / tot as f64;
        }
        Self {
            kmer: self.kmer,
            last: self.last,
            weight,
            probs,
            tot,
            edges,
        }
    }
    fn finalize(&mut self) {
        for i in 0..4 {
            self.probs[i] = self.weight[i] as f64 / self.tot as f64;
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
        self.tot += 1;
        self.weight[i] += 1;
    }
    // return how many occurence there is.
    #[allow(dead_code)]
    fn count(&self) -> u16 {
        self.tot
    }
    // return Score(self.kmer,tip)
    fn calc_score(&self, tip: &[u8], _config: &Config) -> i32 {
        tip.len() as i32 - edlib_sys::global_dist(&self.kmer, tip) as i32
    }
    // return P(idx|self)
    fn to(&self, idx: usize) -> f64 {
        unsafe { *self.probs.get_unchecked(idx) }
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
            b'A' => self.probs[0],
            b'C' => self.probs[1],
            b'G' => self.probs[2],
            b'T' => self.probs[3],
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
        let mut kmer = Kmer::new(b"ACTCGTA");
        assert_eq!(kmer.last, b'A');
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'C', 34);
        kmer.push_edge_with(b'G', 32);
        kmer.push_edge_with(b'A', 12);
        if PSEUDO_COUNT {
            assert_eq!(kmer.weight, [3, 2, 2, 1]);
        } else {
            assert_eq!(kmer.weight, [2, 1, 1, 0]);
        }
        assert_eq!(kmer.edges, [Some(12), Some(34), Some(32), None]);
    }
    #[test]
    #[should_panic]
    fn kmer2() {
        let mut kmer = Kmer::new(b"ACTCGTA");
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'A', 34);
    }
    #[test]
    fn initialize() {
        let test = [
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
        {
            let likelihood1 = model1.forward(&template1, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&template1, &DEFAULT_CONFIG);
            assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
        }
        {
            let likelihood1 = model1.forward(&template2, &DEFAULT_CONFIG);
            let likelihood2 = model2.forward(&template2, &DEFAULT_CONFIG);
            assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
        }
    }
    #[test]
    fn low_coverage_test() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let k = 7;
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
        let model2 = DBGHMM::new(&model2, k);
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
    fn forward_bench(b: &mut Bencher) {
        let bases = b"ACTG";
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let num = 25;
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
                    .map(|e| (e.calc_score(&template[..k], &DEFAULT_CONFIG) as f64).exp())
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
