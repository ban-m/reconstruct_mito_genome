use super::base_table::BASE_TABLE;
// use super::LAMBDA;
use std::fmt;
#[derive(Default, Clone)]
pub struct Base {
    pub base: u8,
    pub edges: Vec<usize>,
    pub weights: Vec<f64>,
    pub base_count: [f64; 4],
    pub weight: f64,
    pub head_weight: f64,
    pub is_tail: bool,
    pub is_head: bool,
}

impl Base {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            edges: vec![],
            weights: vec![],
            weight: 0.,
            head_weight: 0.,
            base_count: [0.; 4],
            is_tail: false,
            is_head: false,
        }
    }
    pub fn weight(&self) -> f64 {
        self.weight
    }
    pub fn head_weight(&self) -> f64 {
        self.head_weight
    }
    pub fn remove_if(&mut self, mapping: &[(usize, bool)]) {
        self.weights = self
            .weights
            .iter()
            .enumerate()
            .filter(|&(idx, _)| mapping[self.edges[idx]].1)
            .map(|(_, &w)| w)
            .collect();
        self.edges = self
            .edges
            .iter()
            .filter_map(|&idx| {
                if mapping[idx].1 {
                    Some(mapping[idx].0)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(self.edges.len(), self.weights.len());
    }
    pub fn remove_lightweight_edges(&mut self) {
        if self.edges.len() > 1 {
            let tot = self.weights.iter().sum::<f64>();
            self.edges = self
                .edges
                .iter()
                .enumerate()
                .filter_map(|(idx, &to)| {
                    if self.weights[idx] / tot > 0.2 {
                        Some(to)
                    } else {
                        None
                    }
                })
                .collect();
            self.weights = self
                .weights
                .iter()
                .filter(|&&w| w / tot > 0.2)
                .copied()
                .collect();
        }
    }
    pub fn finalize(&mut self) {
        // let tot = self.base_count.iter().sum::<f64>();
        // if tot > 0.001 {
        //     self.base_count.iter_mut().for_each(|e| *e /= tot);
        // } else {
        //     self.base_count.iter_mut().for_each(|e| *e = 0.25);
        // }
        let tot = self.weights.iter().sum::<f64>();
        self.weights.iter_mut().for_each(|e| *e /= tot);
    }
    pub fn add(&mut self, b: u8, w: f64, idx: usize) {
        let pos = match self
            .edges
            .iter()
            .enumerate()
            .filter(|&(_, &to)| to == idx)
            .nth(0)
        {
            Some((pos, _)) => pos,
            None => {
                self.edges.push(idx);
                self.weights.push(0.);
                self.edges.len() - 1
            }
        };
        self.weights[pos] += w;
        self.base_count[BASE_TABLE[b as usize]] += w;
    }
    pub fn add_weight(&mut self, w: f64) {
        self.weight += w;
    }
    pub fn add_head_weight(&mut self, w: f64) {
        self.head_weight += w;
    }
    pub fn rename_by(&mut self, map: &[usize]) {
        self.edges.iter_mut().for_each(|e| *e = map[*e]);
    }
    pub fn base(&self) -> u8 {
        self.base
    }
    pub fn to(&self, to: usize) -> f64 {
        *self
            .edges
            .iter()
            .zip(self.weights.iter())
            .filter(|&(&idx, _)| idx == to)
            .nth(0)
            .unwrap()
            .1
    }
    pub fn prob(&self, base: u8, config: &super::Config) -> f64 {
        if self.base == base {
            1. - config.mismatch
        } else {
            config.mismatch / 3.
        }
    }
    #[inline]
    pub fn insertion(&self, _base: u8) -> f64 {
        0.25
        //let q = 0.25;
        // let p = self.base_count[BASE_TABLE[base as usize]];
        //q
        // if self.edge_num() <= 1 {
        //     p * LAMBDA + (1. - LAMBDA) * q
        // } else {
        //     q
        // }
    }
    #[inline]
    pub fn has_edge(&self) -> bool {
        !self.edges.is_empty()
    }
    // #[inline]
    // pub fn edge_num(&self) -> u8 {
    //     self.edges.len() as u8
    // }
    pub fn edges(&self) -> &[usize] {
        &self.edges
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Base\t{}", self.base as char)?;
        let weights: Vec<_> = self.weights.iter().map(|x| format!("{:.3}", x)).collect();
        write!(f, "{}", weights.join("\t"))?;
        for to in self.edges.iter() {
            writeln!(f, "Edge\t{}", to)?;
        }
        writeln!(f, "Is tail\t{}", self.is_tail)?;
        write!(f, "Is head\t{}", self.is_head)
    }
}

impl fmt::Debug for Base {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}\t{}/{}\t{:.2}\t{:.2}",
            self.base as char, self.is_head, self.is_tail, self.weight, self.head_weight,
        )?;
        for (w, to) in self.weights.iter().zip(self.edges.iter()) {
            writeln!(f, "E\t{}\t{:.3}", to, w)?;
        }
        for &b in b"ATGC" {
            let count = self.base_count[BASE_TABLE[b as usize]];
            if count > 0.001 {
                write!(f, "{}:{} ", b as char, count)?;
            }
        }
        Ok(())
    }
}
