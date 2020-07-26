use super::find_union::FindUnion;
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
    cluster: usize,
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
            "{}\t{}\t{}\t[{}]",
            self.cluster,
            self.occ,
            edges.join(","),
            kmer.join(",")
        )
    }
}

impl Node {
    fn new(w: &[(u64, u64)]) -> Self {
        let first = {
            let &(unit, cluster) = w.first().unwrap();
            (unit, cluster)
        };
        let last = {
            let &(unit, cluster) = w.last().unwrap();
            (unit, cluster)
        };
        let kmer: Vec<_> = if first < last {
            w.iter().map(|&(unit, cluster)| (unit, cluster)).collect()
        } else {
            w.iter()
                .rev()
                .map(|&(unit, cluster)| (unit, cluster))
                .collect()
        };
        let (edges, occ, cluster) = (vec![], 0, 0);
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

impl DeBruijnGraph {
    pub fn new(reads: &[Vec<(u64, u64)>], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            for w in read.windows(k + 1) {
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
    pub fn assign_read(&self, read: &[(u64, u64)]) -> Option<usize> {
        let mut count = HashMap::<_, u32>::new();
        for w in read.windows(self.k) {
            let node = Node::new(w);
            if let Some(&idx) = self.indexer.get(&node) {
                *count.entry(self.nodes[idx].cluster).or_default() += 1;
            }
        }
        count.into_iter().max_by_key(|x| x.1).map(|x| x.0)
    }
    pub fn coloring(&mut self) {
        // Coloring node of the de Bruijn graph.
        // As a first try, I just color de Bruijn graph by its connected components.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut current_component = 1;
        debug!("ClusterID\tNumberOfKmer");
        let mut ignored = 0;
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            if count < 10 {
                ignored += 1;
                continue;
            }
            debug!("{}\t{}", current_component, count);
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = current_component;
                }
            }
            current_component += 1;
        }
        debug!("Ignored small component (<10 kmer):{}", ignored);
    }
}
