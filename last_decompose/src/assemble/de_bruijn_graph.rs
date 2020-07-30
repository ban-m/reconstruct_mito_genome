use crate::find_union::FindUnion;
use std::collections::HashMap;
#[derive(Clone)]
pub struct DeBruijnGraph {
    k: usize,
    nodes: Vec<Node>,
    indexer: HashMap<Node, usize>,
}

impl std::fmt::Debug for DeBruijnGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (idx, node) in self.nodes.iter().enumerate() {
            writeln!(f, "{}\t{:?}", idx, node)?;
        }
        write!(f, "K:{}", self.k)
    }
}

#[derive(Clone)]
struct Node {
    occ: usize,
    edges: Vec<Edge>,
    kmer: Vec<(u64, u64)>,
    cluster: Option<usize>,
}

impl std::fmt::Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edges: Vec<_> = self
            .edges
            .iter()
            .map(|e| format!("(->{},{})", e.to, e.weight))
            .collect();
        let kmer: Vec<_> = self
            .kmer
            .iter()
            .map(|(u, c)| format!("{}-{}", u, c))
            .collect();
        write!(
            f,
            "{:?}\t{}\t{}\t[{}]",
            self.cluster,
            self.occ,
            edges.join(","),
            kmer.join(",")
        )
    }
}

impl Node {
    fn new(w: &[super::chunked_read::Node]) -> Self {
        let first = w.first().unwrap().get_tuple();
        let last = w.last().unwrap().get_tuple();
        use super::chunked_read::Node;
        let kmer: Vec<_> = if first < last {
            w.iter().map(Node::get_tuple).collect()
        } else {
            w.iter().rev().map(Node::get_tuple).collect()
        };
        let (edges, occ, cluster) = (vec![], 0, None);
        Self {
            kmer,
            edges,
            occ,
            cluster,
        }
    }
    fn push(&mut self, to: usize) {
        match self.edges.iter_mut().find(|e| e.to == to) {
            Some(x) => {
                x.weight += 1;
            }
            None => self.edges.push(Edge { to, weight: 1 }),
        }
    }
    fn remove_edge(&mut self, to: usize) {
        self.edges.retain(|x| x.to != to);
    }
}

#[derive(Debug, Clone)]
struct Edge {
    to: usize,
    weight: u64,
}

use std::hash::Hasher;
impl std::hash::Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        assert!(!self.kmer.is_empty());
        // Normalize and hashing.
        if self.kmer.first().unwrap() < self.kmer.last().unwrap() {
            for (unit, cluster) in self.kmer.iter() {
                unit.hash(state);
                cluster.hash(state);
            }
        } else {
            for (unit, cluster) in self.kmer.iter().rev() {
                unit.hash(state);
                cluster.hash(state);
            }
        }
    }
}

impl std::cmp::PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        assert!(!self.kmer.is_empty());
        assert!(!other.kmer.is_empty());
        if self.kmer.len() != other.kmer.len() {
            return false;
        }
        let is_self_normed = self.kmer.first().unwrap() < self.kmer.last().unwrap();
        let is_other_normed = other.kmer.first().unwrap() < other.kmer.last().unwrap();
        match (is_self_normed, is_other_normed) {
            (false, false) | (true, true) => self.kmer.iter().zip(&other.kmer).all(|(x, y)| x == y),
            (false, true) | (true, false) => {
                self.kmer.iter().rev().zip(&other.kmer).all(|(x, y)| x == y)
            }
        }
    }
}
impl std::cmp::Eq for Node {}

#[derive(Debug, Clone, Copy)]
struct Bubble {
    // Index of bubble.
    branches: (usize, usize),
    shoot: usize,
    // Root. Where this bubble collapse.
    root: usize,
}

impl DeBruijnGraph {
    pub fn new(reads: &[super::ChunkedRead], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            for w in read.nodes.windows(k + 1) {
                // Calc kmer
                let from = Node::new(&w[..k]);
                let to = Node::new(&w[1..]);
                // Check entry.
                let from = if !indexer.contains_key(&from) {
                    indexer.insert(from.clone(), nodes.len());
                    nodes.push(from);
                    nodes.len() - 1
                } else {
                    *indexer.get(&from).unwrap()
                };
                let to = if !indexer.contains_key(&to) {
                    indexer.insert(to.clone(), nodes.len());
                    nodes.push(to);
                    nodes.len() - 1
                } else {
                    *indexer.get(&to).unwrap()
                };
                nodes[from].occ += 1;
                nodes[to].occ += 1;
                nodes[from].push(to);
                nodes[to].push(from);
            }
        }
        Self { k, nodes, indexer }
    }
    fn calc_thr_edge(&self) -> u64 {
        let counts = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().fold(0, |x, w| x + w.weight))
            .sum::<u64>();
        let len: usize = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        counts / len as u64 / 3
    }
    pub fn clean_up_auto(self) -> Self {
        let thr = self.calc_thr_edge();
        debug!("Removing edges with weight less than {}", thr);
        self.clean_up_2(thr)
    }
    fn clean_up_2(mut self, thr: u64) -> Self {
        // Removing weak and non-solo edge.
        // This is because, solo-edge would be true edge
        // Even if the weight is small.
        // It is not a problem to leave consective false-edge ,
        // as such a cluster would be removed by filtering
        // clusters by their sizes.
        self.nodes
            .iter_mut()
            .filter(|node| node.edges.len() > 2)
            .for_each(|node| {
                // Remove weak edges.
                node.edges.retain(|edge| edge.weight > thr);
            });
        self
    }
    pub fn resolve_bubbles(&mut self, reads: &mut [super::ChunkedRead]) {
        let mut bubble_spec: Vec<_> = self.enumerate_bubbles(reads);
        let mut to_be_removed = vec![false; self.nodes.len()];
        let len = self.nodes.len();
        for index in 0..len {
            if let Some(bubble) = bubble_spec[index] {
                bubble_spec[index] = None;
                // search nearest bubble.
                let pair_pos =
                    match self.search_nearest_bubble(bubble.root, bubble.shoot, &bubble_spec) {
                        Ok(res) => res,
                        Err(is_root) if is_root => {
                            // We cut the bubble arbitrary.
                            let branch = bubble.branches.0;
                            self.nodes[branch].remove_edge(bubble.shoot);
                            self.nodes[bubble.shoot].remove_edge(branch);
                            continue;
                        }
                        Err(_) => continue,
                    };
                // Never panic.
                let pair = bubble_spec[pair_pos].unwrap();
                bubble_spec[index] = None;
                let resolved_pairs = self.resolve_bubble(reads, bubble, pair);
                // First, remove all the edges from branches to the root.
                for bbl in &[pair, bubble] {
                    let ((b0, b1), shoot) = (bbl.branches, bbl.shoot);
                    self.nodes[b0].remove_edge(shoot);
                    self.nodes[b1].remove_edge(shoot);
                    self.nodes[shoot].remove_edge(b0);
                    self.nodes[shoot].remove_edge(b1);
                }
                // Then, connect anchored pairs.
                for (n0, n1) in resolved_pairs {
                    self.nodes[n0].push(n1);
                    self.nodes[n1].push(n0);
                }
                // Lastly, register all the nodes between two shoots as `removed`
                let (mut prev, mut current) = (bubble.root, bubble.shoot);
                while current != pair.root {
                    to_be_removed[current] = true;
                    let next_nodes = &self.nodes[current].edges;
                    assert!(next_nodes.len() == 2);
                    let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                    prev = current;
                    current = next;
                }
                to_be_removed[bubble.root] = true;
                to_be_removed[pair.root] = true;
            }
        }
        // Removing nodes.
        self.remove_nodes(&to_be_removed);
    }
    fn resolve_bubble(
        &self,
        reads: &[super::ChunkedRead],
        bubble0: Bubble,
        bubble1: Bubble,
    ) -> Vec<(usize, usize)> {
        // [00->10, 01->10, 00->11, 01->11];
        let mut connection_counts = [0; 4];
        for read in reads.iter() {
            let (mut b00, mut b01) = (false, false);
            let (mut b10, mut b11) = (false, false);
            for w in read.nodes.windows(self.k) {
                let node = *self.indexer.get(&Node::new(&w)).unwrap();
                b00 = node == bubble0.branches.0;
                b01 = node == bubble0.branches.1;
                b10 = node == bubble1.branches.0;
                b11 = node == bubble1.branches.1;
            }
            assert!(!b00 || !b01);
            assert!(!b10 || !b11);
            if (b00 || b01) && (b10 || b11) {
                let b0 = if b00 { 0 } else { 1 };
                let b1 = if b10 { 0 } else { 1 };
                connection_counts[(b1 << 1) + b0] += 1;
            }
        }
        connection_counts.iter_mut().for_each(|x| {
            if *x < 10 {
                *x = 0;
            }
        });
        // Case0: bubble0.0 -> bubble1.0 and bubble0.1 -> bubble1.1
        let case0 = connection_counts[0] + connection_counts[3];
        // Case1: bubble0.0 -> bubble1.1 and bubble0.0 -> bubble1.0
        let case1 = connection_counts[1] + connection_counts[2];
        let mut pairs = vec![];
        if case0 > case1 {
            if connection_counts[0] > 10 {
                pairs.push((bubble0.branches.0, bubble1.branches.0));
            }
            if connection_counts[3] > 10 {
                pairs.push((bubble0.branches.1, bubble1.branches.1));
            }
        } else {
            if connection_counts[1] > 10 {
                pairs.push((bubble0.branches.1, bubble1.branches.0));
            }
            if connection_counts[2] > 10 {
                pairs.push((bubble0.branches.0, bubble1.branches.1));
            }
        }
        pairs
    }
    fn enumerate_bubbles(&self, reads: &[super::ChunkedRead]) -> Vec<Option<Bubble>> {
        // Enumerate bubble.
        let mut edge_counts: Vec<_> = (0..self.nodes.len()).map(|idx| vec![0; idx]).collect();
        for read in reads.iter() {
            for w in read.nodes.windows(self.k + 1) {
                let from = *self.indexer.get(&Node::new(&w[..self.k])).unwrap();
                let to = *self.indexer.get(&Node::new(&w[1..])).unwrap();
                edge_counts[from.max(to)][from.min(to)] += 1;
            }
        }
        self.nodes
            .iter()
            .enumerate()
            .map(|(shoot, node)| {
                if node.edges.len() != 3 {
                    return None;
                }
                let (to0, to1, to2) = (node.edges[0].to, node.edges[1].to, node.edges[2].to);
                let bet0and1 = edge_counts[to0.max(to1)][to0.min(to1)];
                let bet0and2 = edge_counts[to0.max(to2)][to0.min(to2)];
                let bet1and2 = edge_counts[to1.max(to2)][to1.min(to2)];
                let (branches, root) = if bet0and1 == 0 && bet0and2 > 0 && bet1and2 > 0 {
                    ((to0, to1), to2)
                } else if bet0and1 > 0 && bet0and2 == 0 && bet1and2 > 0 {
                    ((to0, to2), to1)
                } else if bet0and1 > 0 && bet0and2 > 0 && bet1and2 == 0 {
                    ((to1, to2), to1)
                } else {
                    return None;
                };
                Some(Bubble {
                    branches,
                    shoot,
                    root,
                })
            })
            .collect()
    }
    fn remove_nodes(&mut self, to_be_removed: &[bool]) {
        let mut next_index = vec![];
        {
            let mut index = 0;
            for &b in to_be_removed.iter() {
                next_index.push(index);
                index += b as usize;
            }
        }
        let mut index = 0;
        self.nodes.retain(|_| {
            index += 1;
            to_be_removed[index - 1]
        });
        self.nodes.iter_mut().for_each(|n| {
            n.edges.iter_mut().for_each(|x| x.to = next_index[x.to]);
        });
        self.indexer
            .iter_mut()
            .for_each(|(_, x)| *x = next_index[*x]);
    }
    fn search_nearest_bubble(
        &self,
        root: usize,
        shoot: usize,
        bubbles: &[Option<Bubble>],
    ) -> Result<usize, bool> {
        let (mut prev, mut current) = (root, shoot);
        while bubbles[current].is_none() {
            let next_nodes = &self.nodes[current].edges;
            if next_nodes.len() == 1 {
                // Fail to find pair. But this is bubble is terminating.
                return Err(true);
            } else if next_nodes.len() == 2 {
                // Proceed.
                let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                prev = current;
                current = next;
            } else {
                // Fail to find pair. We've reached very complex bubble.
                return Err(false);
            }
        }
        Ok(current)
    }
    pub fn coloring(&mut self) -> usize {
        // Coloring node of the de Bruijn graph.
        // As a first try, I just color de Bruijn graph by its connected components.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut current_component = 0;
        debug!("ClusterID\tNumberOfKmer");
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            if count < 5 {
                continue;
            }
            debug!("{}\t{}", current_component, count);
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = Some(current_component);
                }
            }
            current_component += 1;
        }
        current_component
    }
    pub fn assign_read(&self, read: &super::ChunkedRead) -> Option<u8> {
        let mut count = HashMap::<_, u32>::new();
        for w in read.nodes.windows(self.k) {
            let node = Node::new(w);
            if let Some(&idx) = self.indexer.get(&node) {
                if let Some(cl) = self.nodes[idx].cluster {
                    *count.entry(cl).or_default() += 1;
                }
            }
        }
        // If this is Some(res), then, x.1 never be 0.
        count.into_iter().max_by_key(|x| x.1).map(|x| x.0 as u8)
    }
}
