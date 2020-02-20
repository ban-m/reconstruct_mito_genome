#![feature(test)]
extern crate log;
extern crate packed_simd;
extern crate rand;
extern crate rand_xoshiro;
extern crate test;
use packed_simd::f64x4 as f64s;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
mod config;
pub use config::*;
mod base;
mod base_table;
use base::Base;
const LAMBDA: f64 = 0.3;
const SMALL: f64 = 0.000_000_001;

// Edit operation
#[derive(Debug, Clone, Copy)]
enum EditOp {
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

use std::fmt;
impl fmt::Debug for POA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let result: Vec<_> = self.nodes.iter().map(|e| format!("{:?}", e)).collect();
        write!(f, "{}", result.join("\n"))
    }
}
impl fmt::Display for POA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let edges = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().count())
            .sum::<usize>();
        write!(f, "{}\t{}", self.nodes.len(), edges)
    }
}

impl PartialOrderAlignment {
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
        Self { weight, nodes }
    }
    fn score(x: u8, y: u8) -> f64 {
        if x == y {
            1.
        } else {
            -1.
        }
    }
    pub fn add(mut self, seq: &[u8], w: f64) -> Self {
        if self.weight < SMALL {
            return Self::new(seq, w);
        }
        // Alignment
        let mut dp: Vec<Vec<_>> = vec![vec![0.; seq.len() + 1]; self.nodes.len() + 1];
        let mut traceback = vec![vec![EditOp::Stop; seq.len() + 1]; self.nodes.len() + 1];
        let edges = self.edges();
        assert!(self.nodes.iter().any(|e| e.is_head));
        let (ins, del) = (-1., -1.);
        for (idx, (dist, parent)) in self.min_dist().into_iter().enumerate() {
            dp[idx + 1][0] = del * dist as f64;
            traceback[idx + 1][0] = EditOp::Deletion(parent);
        }
        for j in 0..seq.len() {
            dp[0][j + 1] = ins * (j + 1) as f64;
            traceback[0][j + 1] = EditOp::Insertion(j);
        }
        for (i, n) in self.nodes.iter().enumerate() {
            for (j, &b) in seq.iter().enumerate() {
                // Match state. j-1 should be valid, because if j == 0, the node should be
                // a root node with no parents.
                let ms = Self::score(n.base(), b);
                let (m_argmax, m_max) = edges[i]
                    .iter()
                    .map(|&from| (from + 1, dp[from + 1][j] + ms))
                    .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
                    .unwrap_or((0, ms));
                // Deletion from the reference(graph)
                let (d_argmax, d_max) = edges[i]
                    .iter()
                    .map(|&from| (from + 1, dp[from + 1][j + 1] + del))
                    .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
                    .unwrap_or((0, -100_000.));
                // Insertion into the reference(graph)
                let (i_argmax, i_max) = (i + 1, dp[i + 1][j] + ins);
                // Select one of the optimal operations.
                let (argmax, max) = if d_max < m_max && i_max < m_max {
                    (EditOp::Match(m_argmax), m_max)
                } else if m_max < d_max && i_max < d_max {
                    (EditOp::Deletion(d_argmax), d_max)
                } else {
                    (EditOp::Insertion(i_argmax), i_max)
                };
                // Fill.
                dp[i + 1][j + 1] = max;
                traceback[i + 1][j + 1] = argmax;
            }
        }
        // Traceback
        // g_pos = position on the graph, q_pos = position on the query
        let mut q_pos = seq.len();
        let (mut g_pos, _) = self
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(idx, n)| {
                if n.is_tail {
                    Some((idx + 1, dp[idx + 1][q_pos]))
                } else {
                    None
                }
            })
            .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap())
            .unwrap();
        let mut previous = None;
        while q_pos > 0 {
            match traceback[g_pos][q_pos] {
                EditOp::Match(from) => {
                    let current_pos = if self.nodes[g_pos - 1].base() == seq[q_pos - 1] {
                        g_pos - 1
                    } else {
                        self.nodes.push(Base::new(seq[q_pos - 1]));
                        self.nodes.len() - 1
                    };
                    if let Some(p) = previous {
                        let base = self.nodes[p as usize].base();
                        self.nodes[current_pos].add(base, w, p);
                    };
                    previous = Some(current_pos);
                    g_pos = from;
                    q_pos -= 1;
                }
                EditOp::Deletion(from) => {
                    g_pos = from;
                }
                EditOp::Insertion(_) => {
                    self.nodes.push(Base::new(seq[q_pos - 1]));
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
        self.select_head_tail();
        self.topological_sort()
    }
    // Return the minimum distance from root node to each node.
    // Done by BFS.
    fn min_dist(&self) -> Vec<(usize, usize)> {
        let mut dist = vec![(0, 0); self.nodes.len()];
        let mut queue: Vec<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| if e.is_head { Some(idx) } else { None })
            .collect();
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
    fn select_head_tail(&mut self) {
        let mut is_head = vec![true; self.nodes.len()];
        for n in self.nodes.iter() {
            for &to in n.edges.iter() {
                is_head[to] = true;
            }
        }
        self.nodes.iter_mut().zip(is_head).for_each(|(n, is_head)| {
            n.is_tail = n.edges.is_empty();
            n.is_head = is_head;
        });
    }
    fn topological_sort(mut self) -> Self {
        let mapping = Self::topological_sort_inner(&self.nodes);
        let mut result = vec![None; self.nodes.len()];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            node.rename_by(&mapping);
            result[mapping[idx]] = Some(node);
        }
        assert!(result.iter().all(|e| e.is_some()));
        assert!(self.nodes.is_empty());
        self.nodes.extend(result.into_iter().filter_map(|e| e));
        self.nodes
            .iter_mut()
            .for_each(|n| n.is_tail = n.edges.is_empty());
        self
    }
    // Return topological sorted order. The graph never contains cycles.
    fn topological_sort_inner(nodes: &[Base]) -> Vec<usize> {
        let mut edges = vec![vec![]; nodes.len()];
        for (idx, node) in nodes.iter().enumerate() {
            for &to in &node.edges {
                edges[idx].push(to);
            }
        }
        Self::topological_sort_dfs(&edges)
    }
    // Topological sorting.
    fn topological_sort_dfs(edges: &[Vec<usize>]) -> Vec<usize> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let len = edges.len();
        let mut dfs_flag = vec![0; len];
        let mut dfs_stack = vec![];
        let mut mapping = vec![];
        for i in 0..len {
            if dfs_flag[i] != 0 {
                continue;
            }
            assert!(dfs_stack.is_empty());
            dfs_stack.push(i);
            'dfs: while !dfs_stack.is_empty() {
                let node = *dfs_stack.last().unwrap();
                if dfs_flag[node] == 0 {
                    // preorder
                    dfs_flag[node] = 1;
                }
                for &to in &edges[node] {
                    if dfs_flag[to] == 0 {
                        dfs_stack.push(to);
                        continue 'dfs;
                    }
                }
                // No-op
                let last = dfs_stack.pop().unwrap();
                mapping.push(last);
                // Deactivate
                dfs_flag[last] = 2;
            }
        }
        mapping.reverse();
        assert!(dfs_stack.is_empty());
        // Use DFS stack as a temporary buffer.
        std::mem::swap(&mut mapping, &mut dfs_stack);
        // Filling zero so that the boundary would be violated.
        mapping.extend(std::iter::repeat(0).take(len));
        for (order, &v) in dfs_stack.iter().enumerate() {
            mapping[v] = order;
        }
        mapping
    }
    pub fn edges(&self) -> Vec<Vec<usize>> {
        let mut edges = vec![vec![]; self.nodes.len()];
        for (from, n) in self.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[to].push(from);
            }
        }
        edges
    }
    pub fn clean_up(mut self) -> Self {
        self.nodes.iter_mut().for_each(|e| e.finalize());
        self
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
    /// Return the total weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    fn update(
        &self,
        updates: &mut [f64],
        prev: &[f64],
        base: u8,
        config: &Config,
        edges: &[Vec<usize>],
    ) -> (f64, f64) {
        assert!((1. - prev.iter().sum::<f64>()).abs() < SMALL);
        // Alignemnt:[mat,ins,del, mat,ins,del, mat,....,del]
        for (dist_idx, (dist, froms)) in self.nodes.iter().zip(edges.iter()).enumerate() {
            let node = 3 * dist_idx;
            let (match_state, insertion_state) = froms
                .iter()
                .map(|&src| {
                    let src_node = &self.nodes[src];
                    let src = src * 3;
                    let trans = src_node.to(dist_idx);
                    let f_dist = prev[src] * config.p_match
                        + prev[src + 1] * (1. - config.p_extend_ins)
                        + prev[src + 2] * (1. - config.p_extend_del - config.p_del_to_ins);
                    let m = f_dist * trans * dist.prob(base, config);
                    let i = prev[node + 2] * config.p_del_to_ins * trans;
                    (m, i)
                })
                .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
            updates[node] = match_state;
            updates[node + 1] = insertion_state;
            updates[node + 1] += if dist.has_edge() {
                prev[node] * config.p_ins + prev[node + 1] * config.p_extend_ins
            } else {
                prev[node..=node + 2].iter().sum::<f64>()
            };
            updates[node + 1] *= dist.insertion(base);
        }
        let d = Self::sum(updates).recip();
        Self::mul(updates, d);
        for (dist_idx, src_nodes) in edges.iter().enumerate() {
            updates[3 * dist_idx + 2] = src_nodes
                .iter()
                .map(|&src| {
                    let trans = updates[3 * src] * config.p_del
                        + updates[3 * src + 2] * config.p_extend_del;
                    trans * self.nodes[src].to(dist_idx)
                })
                .sum::<f64>();
        }
        let c = Self::sum(&updates).recip();
        Self::mul(updates, c);
        (c, d)
    }
    #[cfg(target_feature = "sse")]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        // Alignemnts: [mat, ins, del,  mat, ins, del,  ....]
        let heads = self.nodes.iter().filter(|n| n.is_head).count() as f64;
        let mut prev: Vec<_> = self
            .nodes
            .iter()
            .flat_map(|n| {
                if n.is_head {
                    vec![1. / heads, 0., 0.]
                } else {
                    vec![0.; 3]
                }
            })
            .collect();
        if prev.len() % 4 != 0 {
            prev.extend(std::iter::repeat(0.).take(4 - prev.len() % 4));
        }
        let mut lk = 0.;
        let mut updated = vec![0.; prev.len()];
        let edges = self.edges();
        for (idx, &base) in obs.iter().enumerate() {
            updated.iter_mut().for_each(|e| *e = 0.);
            let (c, d) = self.update(&mut updated, &prev, base, config, &edges);
            assert!(c.ln() + d.ln() > 0., "{},{}", c, d);
            lk -= if idx + 1 < obs.len() {
                c.ln() + d.ln()
            } else {
                d.ln()
            };
            eprintln!("{}\t{}", idx, lk);
            std::mem::swap(&mut prev, &mut updated);
        }
        lk
    }
}

pub fn generate(seqs: &[&[u8]], ws: &[f64]) -> POA {
    let seed = ws.iter().sum::<f64>().floor() as u64 + seqs.len() as u64;
    let choises: Vec<_> = (0..seqs.len()).collect();
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let picked = *choises.choose_weighted(&mut rng, |&i| ws[i]).unwrap();
    let (seed, seed_weight) = (&seqs[picked], ws[picked]);
    seqs.iter()
        .zip(ws.iter())
        .filter(|&(_, &w)| w > 0.001)
        .fold(POA::new(seed, seed_weight / 2.), |x, (y, &w)| x.add(y, w))
        .clean_up()
}

pub fn generate_uniform(seqs: &[&[u8]]) -> POA {
    let ws = vec![1.; seqs.len()];
    generate(seqs, &ws)
}

pub fn generate_vec(seqs: &[Vec<u8>]) -> POA {
    let ws = vec![1.; seqs.len()];
    let seqs: Vec<_> = seqs.iter().map(|e| e.as_slice()).collect();
    generate(&seqs, &ws)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
    #[test]
    fn init() {
        let _poa = POA::default();
        let seed = b"ACGTAGCTGATCGTAC";
        let _poa = POA::new(seed, 1.);
        let poa = POA::new(seed, 0.4);
        assert!(poa.nodes.iter().any(|n| n.is_head));
    }
    #[test]
    fn add() {
        let seed = b"ACGTAGCTGATCGTAC";
        POA::new(seed, 1.).add(seed, 1.);
    }
    #[test]
    fn create() {
        let seed = b"ACGTAGCTGATCGTAC";
        let res = POA::new(seed, 1.).add(seed, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len());
    }
    #[test]
    fn bubble() {
        let seed = b"ACGTAGCTGATCGTAC";
        let seq = b"ACGTAGCTGATCGGAC";
        let res = POA::new(seed, 1.).add(seq, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 1);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATCGTCC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 1);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATTTCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 2);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATTTCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 2);
        let seed = b"ACGTAGCTGAGTAC";
        let seq2 = b"ACGAGCTGATCTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 2);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len());
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATCGTACG";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 1);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"CCGTAGCTGATCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 1);
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"ACGTAGCTGATCGTAG";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 1);
        let seq2 = b"ACGTAGCTGATCGTAC";
        let seed = b"AACGTAGCTGATCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len());
        let seed = b"ACGTAGCTGATCGTAC";
        let seq2 = b"AAACGTAGCTGATCGTAC";
        let res = POA::new(seed, 1.).add(seq2, 1.).clean_up();
        assert_eq!(res.nodes.len(), seed.len() + 2);
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
        test.iter()
            .fold(POA::default(), |x, y| x.add(y, 1.))
            .clean_up();
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
        let m = generate_vec(&test);
        eprintln!("{}", m);
        let lk = m.forward(b"CACACAGCAGTCAGTGCA", &DEFAULT_CONFIG);
        eprintln!("{}", lk);
        assert!(lk < 0.);
    }
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
    fn connectivity_check() {
        let bases = b"ACTG";
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
        let template: Vec<_> = (0..30)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..20)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let mut m = POA::default();
        for seq in model1 {
            m = m.add(&seq, 1.);
            assert!(is_connected(&m));
        }
    }
    fn is_connected(m: &POA) -> bool {
        let mut edges = vec![vec![]; m.nodes.len()];
        for (from, n) in m.nodes.iter().enumerate() {
            for &to in n.edges.iter() {
                edges[from].push(to);
                edges[to].push(from);
            }
        }
        let mut dfs = vec![false; edges.len()];
        let mut stack = vec![0];
        'dfs: while !stack.is_empty() {
            let n = *stack.last().unwrap();
            if !dfs[n] {
                dfs[n] = true;
            }
            for &to in edges[n].iter() {
                if !dfs[to] {
                    stack.push(to);
                    continue 'dfs;
                }
            }
            stack.pop();
        }
        dfs.iter().all(|&e| e)
    }
    #[test]
    fn forward_check() {
        let bases = b"ACTG";

        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
        let template: Vec<_> = (0..30)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..20)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let m = generate_vec(&model1);
        let lk = m.forward(&template, &DEFAULT_CONFIG);
        eprintln!("{}", m);
        assert!(lk < 0., "{}", lk)
    }
    #[test]
    fn random_check() {
        let bases = b"ACTG";
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
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
        let model1 = generate_vec(&model1);
        let model2 = generate_vec(&model2);
        let likelihood1 = model1.forward(&template, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&template, &DEFAULT_CONFIG);
        assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
    }
    #[test]
    fn hard_test() {
        let bases = b"ACTG";
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
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
        let model1 = generate_vec(&model1);
        let model2 = generate_vec(&model2);
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
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
        let len = 150;
        let template1: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let p = Profile {
            sub: 0.03,
            ins: 0.03,
            del: 0.03,
        };
        let cov = 30;
        for _ in 0..5 {
            let template2 = introduce_randomness(&template1, &mut rng, &p);
            let model1: Vec<Vec<_>> = (0..cov)
                .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
                .collect();
            let model2: Vec<Vec<_>> = (0..cov)
                .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
                .collect();
            let dataset: Vec<_> = model1
                .iter()
                .map(|e| e.as_slice())
                .chain(model2.iter().map(|e| e.as_slice()))
                .collect();
            let weight1 = vec![vec![0.8; cov], vec![0.2; cov]].concat();
            let weight2 = vec![vec![0.2; cov], vec![0.8; cov]].concat();
            let model1 = generate(&dataset, &weight1);
            let model2 = generate(&dataset, &weight2);
            let num = 50;
            let correct = (0..num)
                .filter(|_| {
                    let q = introduce_randomness(&template1, &mut rng, &PROFILE);
                    let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                    let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                    eprintln!("1\t{:.3}\t{:.3}", lk1, lk2);
                    lk1 > lk2
                })
                .count();
            assert!(correct >= num * 4 / 5, "{}", correct);
            let correct = (0..num)
                .filter(|_| {
                    let q = introduce_randomness(&template2, &mut rng, &PROFILE);
                    let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                    let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                    eprintln!("2\t{:.3}\t{:.3}", lk1, lk2);
                    lk1 < lk2
                })
                .count();
            assert!(correct >= num * 4 / 5, "{}", correct);
        }
    }
    #[test]
    fn abundance_test_prior() {
        let bases = b"ACTG";
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_219);
        let p = Profile {
            sub: 0.03,
            ins: 0.03,
            del: 0.03,
        };
        let len = 150;
        let cov = 20;
        let ratio = 5;
        let errors = PROFILE;
        for _ in 0..3 {
            let template1: Vec<_> = (0..len)
                .filter_map(|_| bases.choose(&mut rng))
                .copied()
                .collect();
            let template2 = introduce_randomness(&template1, &mut rng, &p);
            let data1: Vec<Vec<_>> = (0..(ratio + 1) * cov)
                .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
                .collect();
            let mut data2: Vec<Vec<_>> = (0..ratio * cov)
                .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
                .collect();
            data2.extend((0..cov).map(|_| introduce_randomness(&template2, &mut rng, &PROFILE)));
            let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
            let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
            let total = (ratio + 1) * cov;
            let weight = vec![1.; total];
            let model1 = generate(&data1, &weight);
            let model2 = generate(&data2, &weight);
            let num = 50;
            let correct = (0..num)
                .filter(|_| {
                    let q = introduce_randomness(&template1, &mut rng, &errors);
                    let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                    let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                    eprintln!("1\t{:.3}\t{:.3}", lk1, lk2);
                    lk1 > lk2
                })
                .count();
            eprintln!("1:{}", correct);
            assert!(correct >= num * 6 / 10, "1:{}", correct);
            let correct = (0..num)
                .filter(|_| {
                    let q = introduce_randomness(&template2, &mut rng, &errors);
                    let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                    let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                    eprintln!("2\t{:.3}\t{:.3}", lk1, lk2);
                    lk1 < lk2
                })
                .count();
            eprintln!("2:{}", correct);
            assert!(correct >= num * 6 / 10, "2:{}", correct);
        }
    }
    #[test]
    fn single_error_test() {
        let bases = b"ACTG";
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
        let coverage = 200;
        let start = 20;
        let step = 4;
        let len = 150;
        let results: Vec<_> = (start..coverage)
            .step_by(step)
            .map(|cov| {
                let template1: Vec<_> = (0..len)
                    .filter_map(|_| bases.choose(&mut rng))
                    .copied()
                    .collect();
                let template2 = introduce_errors(&template1, &mut rng, 1, 0, 0);
                let sub = check(&template1, &template2, &mut rng, cov);
                let template2 = introduce_errors(&template1, &mut rng, 0, 1, 0);
                let del = check(&template1, &template2, &mut rng, cov);
                let template2 = introduce_errors(&template1, &mut rng, 0, 0, 1);
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
        eprintln!("Tot:{}", (start..coverage).step_by(step).count() * 100);
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
        let m1 = generate(&seqs, &weight1);
        let m2 = generate(&seqs, &weight2);
        eprintln!("{}\t{}\t{}", cov, m1, m2);
        let correct = (0..100)
            .filter(|e| {
                if e % 2 == 0 {
                    let q = introduce_randomness(&t1, rng, &PROFILE);
                    m1.forward(&q, &DEFAULT_CONFIG) > m2.forward(&q, &DEFAULT_CONFIG)
                } else {
                    let q = introduce_randomness(&t2, rng, &PROFILE);
                    m1.forward(&q, &DEFAULT_CONFIG) < m2.forward(&q, &DEFAULT_CONFIG)
                }
            })
            .count();
        correct
    }
    #[test]
    fn low_coverage_test() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(121212);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
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
        let model1 = generate_vec(&model1);
        eprintln!("Model1:{}", model1);
        let model2 = generate_vec(&model2);
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
    #[test]
    fn low_coverage_weighted_test() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(121212);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..10)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
        let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
        let dataset: Vec<_> = model1
            .iter()
            .chain(model2.iter())
            .map(|e| e.as_slice())
            .collect();
        let weight1 = vec![vec![1.; 10], vec![0.; 10]].concat();
        let weight2 = vec![vec![0.; 10], vec![1.; 10]].concat();
        let model1 = generate(&dataset, &weight1);
        eprintln!("Model1:{}", model1);
        let model2 = generate(&dataset, &weight2);
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

    #[test]
    fn high_coverage_test() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(121212332);
        let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
        let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
        let model1: Vec<Vec<_>> = (0..200)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..200)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
        let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
        {
            let model1 = generate_vec(&model1);
            let model2 = generate_vec(&model2);
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
        let model1 = generate(&dataset, &weight1);
        eprintln!("Model1:{}", model1);
        let model2 = generate(&dataset, &weight2);
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
    #[derive(Debug, Clone, Copy)]
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
    pub fn introduce_errors<T: rand::Rng>(
        seq: &[u8],
        rng: &mut T,
        sub: usize,
        del: usize,
        ins: usize,
    ) -> Vec<u8> {
        // Alignment operations.
        let mut operations = vec![
            vec![Op::Match; seq.len() - sub - del],
            vec![Op::MisMatch; sub],
            vec![Op::Del; del],
            vec![Op::In; ins],
        ]
        .concat();
        operations.shuffle(rng);
        let mut res = vec![];
        let mut remainings: Vec<_> = seq.iter().copied().rev().collect();
        for op in operations {
            match op {
                Op::Match => res.push(remainings.pop().unwrap()),
                Op::MisMatch => res.push(choose_base(rng, remainings.pop().unwrap())),
                Op::In => res.push(random_base(rng)),
                Op::Del => {
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
