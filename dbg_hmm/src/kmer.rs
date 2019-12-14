use super::find_union;
use super::Config;
use super::PSEUDO_COUNT;
#[derive(Clone, Default)]
pub struct Kmer {
    pub kmer: Vec<u8>,
    last: u8,
    // parameters. First four for weight,
    // last four for transition.
    weight: [f64; 8],
    // weight: [f64; 4],
    // transition: [f64; 4],
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
}

impl std::fmt::Debug for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Kmer:{}", String::from_utf8_lossy(&self.kmer))?;
        writeln!(f, "KmerWeight:{:.3}", self.kmer_weight)?;
        writeln!(f, "Last:{}", self.last as char)?;
        // writeln!(f, "Weight:{:?}", &self.weight[..4])?;
        // writeln!(f, "Transition:{:?}", &self.weight[4..])?;
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
    pub fn new(x: &[u8], w: f64) -> Self {
        let kmer = x.to_vec();
        let kmer_weight = w;
        let last = *kmer.last().unwrap();
        // Prior
        let mut weight = [0.; 8];
        if PSEUDO_COUNT {
            for i in 0..4 {
                weight[i] = 1.;
            }
        }
        // let weight = if PSEUDO_COUNT { [1.; 4] } else { [0.; 4] };
        // let transition = [0f64; 4];
        let tot = 0.;
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
            //  transition,
            is_tail,
            is_head,
        }
    }
    pub fn finalize(&mut self) {
        let tot_for_weight = if PSEUDO_COUNT {
            self.tot + 4.
        } else {
            self.tot
        };
        if self.tot > 0. {
            for i in 0..4 {
                self.weight[i] /= tot_for_weight;
                //self.transition[i] /= self.tot;
                self.weight[i + 4] /= self.tot;
            }
        } else {
            for i in 0..4 {
                self.weight[i] = 0.25;
            }
        }
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
                    // self.tot -= self.transition[i];
                    // self.weight[i] -= self.transition[i];
                    // self.transition[i] = 0.;
                    self.tot -= self.weight[i + 4];
                    self.weight[i] -= self.weight[i + 4];
                    self.weight[i + 4] = 0.;
                }
            }
        }
    }
    // Remove all the edges to unsuppoeted nodes.
    pub fn remove_if_not_supported(&mut self, is_supported: &[u8]) {
        for i in 0..4 {
            if let Some(res) = self.edges[i] {
                if is_supported[res] != 1 {
                    self.edges[i] = None;
                    // self.tot -= self.transition[i];
                    // self.weight[i] -= self.transition[i];
                    // self.transition[i] = 0.;
                    self.tot -= self.weight[i + 4];
                    self.weight[i] -= self.weight[i + 4];
                    self.weight[i + 4] = 0.;
                }
            }
        }
    }
    pub fn push_edge_with(&mut self, base: u8, to: usize) {
        self.push_edge_with_weight(base, to, 1.);
    }
    pub fn push_edge_with_weight(&mut self, base: u8, to: usize, w: f64) {
        use super::base_table::BASE_TABLE;
        let i = BASE_TABLE[base as usize];
        self.edges[i] = Some(to);
        self.tot += w;
        self.weight[i] += w;
        //self.transition[i] += w;
        self.weight[i+4] += w;
    }
    // return [mat, ins, del, mism]
    pub fn calc_score(&self, tip: &[u8]) -> [u8; 4] {
        let aln = edlib_sys::global(&self.kmer, tip);
        aln.into_iter().fold([0; 4], |mut xs, x| {
            xs[x as usize] += 1;
            xs
        })
    }
    // return P(idx|self)
    #[inline]
    pub fn to(&self, idx: usize) -> f64 {
        //self.transition[idx]
        self.weight[idx + 4]
    }
    #[inline]
    pub fn prob(&self, base: u8, config: &Config) -> f64 {
        if self.last != base {
            config.mismatch / 3.
        } else {
            1. - config.mismatch
        }
    }
    // return P_I(base|self)
    pub fn insertion(&self, base: u8) -> f64 {
        use super::base_table::BASE_TABLE;
        self.weight[BASE_TABLE[base as usize]]
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
        if PSEUDO_COUNT {
            assert_eq!(&kmer.weight[..4], &[1. + 1. + 1., 1. + 1., 1. + 1., 1.]);
        } else {
            assert_eq!(&kmer.weight[..4], &[1. + 1., 1., 1., 0.]);
        }
        assert_eq!(kmer.edges, [Some(12), Some(34), Some(32), None]);
    }
    #[test]
    fn kmer2() {
        let mut kmer = Kmer::new(b"ACTCGTA", 1.);
        kmer.push_edge_with(b'A', 12);
        kmer.push_edge_with(b'A', 34);
    }
}
