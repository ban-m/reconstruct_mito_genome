use super::base_table::BASE_TABLE;
use super::find_union;
use super::Config;
use super::PSEUDO_COUNT;
const LAMBDA: f64 = 0.10;
#[derive(Clone, Default)]
pub struct Kmer {
    pub kmer: Vec<u8>,
    last: u8,
    // Transition counts.
    pub weight: [f64; 4],
    pub base_count: [f64; 4],
    /// Total number of *Outdegree*
    pub tot: f64,
    /// The location to the edges with the label of A,C,G,and T.
    /// If there is no edges, None
    pub edges: [Option<usize>; 4],
    /// Weight of this kmer.
    pub kmer_weight: f64,
    // Whether this is the end of unit.
    pub is_tail: bool,
    pub is_head: bool,
    pub edge_num: u8,
}

impl std::fmt::Debug for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Kmer:{}", String::from_utf8_lossy(&self.kmer))?;
        writeln!(f, "KmerWeight:{:.3}", self.kmer_weight)?;
        writeln!(f, "Last:{}", self.last as char)?;
        writeln!(f, "Transition:{:?}", &self.weight)?;
        writeln!(f, "BaseCount{:?}", &self.base_count)?;
        writeln!(f, "tot:{}", self.tot)?;
        writeln!(f, "is_tail:{}", self.is_tail)?;
        writeln!(f, "is_head:{}", self.is_head)?;
        write!(f, "edge_num:{}", self.edge_num)?;
        let mut res = String::new();
        for (i, to) in self
            .edges
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| e.map(|to| (idx, to)))
        {
            res.push_str(&format!("\n->{}({})", to, i));
        }
        write!(f, "{}", res)
    }
}

impl Kmer {
    pub fn new(x: &[u8], w: f64) -> Self {
        let kmer = x.to_vec();
        let kmer_weight = w;
        let last = *kmer.last().unwrap();
        // Prior
        let weight = [0.; 4];
        let base_count = [PSEUDO_COUNT; 4];
        let tot = 0.;
        let edges = [None; 4];
        let is_tail = false;
        let is_head = false;
        let edge_num = 0;
        Self {
            kmer,
            kmer_weight,
            last,
            weight,
            base_count,
            tot,
            edges,
            is_tail,
            is_head,
            edge_num,
        }
    }
    pub fn finalize(&mut self) {
        if self.tot > 0.0001 {
            for i in 0..4 {
                self.weight[i] /= self.tot;
            }
        }
        let tot = self.base_count.iter().sum::<f64>();
        self.edge_num = self.edges.iter().filter(|e| e.is_some()).count() as u8;
        if tot > 0.01 {
            self.base_count.iter_mut().for_each(|e| *e /= tot);
        } else {
            self.base_count = [0.25; 4];
        }
        assert!((1. - self.base_count.iter().sum::<f64>()).abs() < 0.001);
    }
    // If which[i] = true, the weight of these edges would be normalized.
    pub fn finalize_global(&mut self, which: &[bool; 4]) {
        assert!(self.tot > 0.0001);
        let tot = self
            .weight
            .iter()
            .zip(which.iter())
            .filter(|&(_, &b)| b)
            .map(|(w, _)| w)
            .sum::<f64>();
        for (i, &b) in which.iter().enumerate() {
            if b {
                self.weight[i] /= tot;
            } else {
                self.weight[i] = if self.weight[i] > 0.0001 { 1. } else { 0. };
            }
        }
        let tot = self.base_count.iter().sum::<f64>();
        self.edge_num = self.edges.iter().filter(|e| e.is_some()).count() as u8;
        if tot > 0.01 {
            self.base_count.iter_mut().for_each(|e| *e /= tot);
        } else {
            self.base_count = [0.25; 4];
        }
        assert!((1. - self.base_count.iter().sum::<f64>()).abs() < 0.001);
    }
    // renaming all the edges by `map`
    pub fn rename_by(&mut self, map: &[usize]) {
        for edge in self.edges.iter_mut() {
            if let Some(res) = edge.as_mut() {
                *res = map[*res];
            }
        }
    }
    // Remove all the edges to nodes not in the maximum group.
    pub fn remove_if_not(&mut self, fu: &mut find_union::FindUnion, mg: usize) {
        for i in 0..4 {
            if let Some(res) = self.edges[i] {
                if fu.find(res).unwrap() != mg {
                    self.edges[i] = None;
                    self.tot -= self.weight[i];
                    self.weight[i] = 0.;
                }
            }
        }
    }
    /// Remove the i-th edge.
    pub fn remove(&mut self, i: usize) {
        self.edges[i] = None;
        self.tot -= self.weight[i];
        self.weight[i] = 0.;
    }
    // Remove all the edges to unsuppoeted nodes.
    pub fn remove_if_not_supported(&mut self, is_supported: &[u8]) {
        for i in 0..4 {
            if let Some(res) = self.edges[i] {
                if is_supported[res] != 1 {
                    self.edges[i] = None;
                    self.tot -= self.weight[i];
                    self.weight[i] = 0.;
                }
            }
        }
    }
    pub fn push_edge_with(&mut self, base: u8, to: usize) {
        self.push_edge_with_weight(base, to, 1.);
    }
    pub fn push_edge_with_weight(&mut self, base: u8, to: usize, w: f64) {
        let i = BASE_TABLE[base as usize];
        self.edges[i] = Some(to);
        self.tot += w;
        self.weight[i] += w;
        self.base_count[i] += w;
    }
    // return P(idx|self)
    #[inline]
    pub fn to(&self, idx: usize) -> f64 {
        self.weight[idx]
    }
    // return P(base|self), observation probability.
    #[inline]
    pub fn prob(&self, base: u8, config: &Config) -> f64 {
        let p = self.base_count[BASE_TABLE[base as usize]];
        let q = if self.last == base {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        };
        p * LAMBDA + (1. - LAMBDA) * q
    }
    #[inline]
    pub fn prob_with(&self, base: u8, config: &Config, from: &Kmer) -> f64 {
        let q = if self.last == base {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        };
        let p = from.base_count[BASE_TABLE[base as usize]];
        if from.edge_num <= 1 {
            p * LAMBDA + (1. - LAMBDA) * q
        } else {
            q
        }
    }
    // return P_I(base|self)
    #[inline]
    pub fn insertion(&self, base: u8) -> f64 {
        let q = 0.25;
        let p = self.base_count[BASE_TABLE[base as usize]];
        if self.edge_num <= 1 {
            p * LAMBDA + (1. - LAMBDA) * q
        } else {
            q
        }
    }
    #[inline]
    pub fn last(&self) -> u8 {
        self.last
    }
    #[inline]
    pub fn has_edge(&self) -> bool {
        self.edge_num != 0
    }
    #[inline]
    pub fn edge_num(&self) -> u8 {
        self.edge_num
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn kmer() {
        let mut kmer = Kmer::new(b"ACTCGTA", 0.);
        assert_eq!(kmer.last, b'A');
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'C', 34);
        kmer.push_edge_with(b'G', 32);
        kmer.push_edge_with(b'A', 12);
        assert_eq!(&kmer.weight, &[1. + 1., 1., 1., 0.]);
        assert_eq!(kmer.edges, [Some(12), Some(34), Some(32), None]);
    }
    #[test]
    fn kmer2() {
        let mut kmer = Kmer::new(b"ACTCGTA", 1.);
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'A', 34);
    }
}
