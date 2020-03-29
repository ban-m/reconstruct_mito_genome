use crate::Config;
use crate::PartialOrderAlignment;
use crate::DEFAULT_LK;
use crate::SMALL;
use packed_simd::f64x4 as f64s;
impl PartialOrderAlignment {
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
    fn update_row(
        &self,
        updates: &mut [f64],
        prev: &[f64],
        base: u8,
        config: &Config,
        edges: &[Vec<(usize, f64)>],
        head_base_freq: &[f64; 4],
    ) -> (f64, f64) {
        use super::base_table::BASE_TABLE;
        assert!((1. - prev.iter().sum::<f64>()).abs() < SMALL);
        let del_to_match = 1. - config.p_extend_del - config.p_del_to_ins;
        for (dist_idx, froms) in edges.iter().enumerate() {
            let match_transition = froms
                .iter()
                .map(|&(src, weight)| {
                    let src_pos = 3 * src;
                    (prev[src_pos] * config.p_match
                        + prev[src_pos + 1] * (1. - config.p_extend_ins)
                        + prev[src_pos + 2] * del_to_match)
                        * weight
                })
                .sum::<f64>();
            use std::cmp::Ordering;
            let match_observation = match dist_idx.cmp(&self.nodes.len()) {
                Ordering::Less => self.nodes[dist_idx].prob(base, config),
                Ordering::Equal => head_base_freq[BASE_TABLE[base as usize]],
                Ordering::Greater => {
                    assert!(match_transition.abs() < 0.00001);
                    1.
                }
            };
            let node = 3 * dist_idx;
            updates[node] = match_transition * match_observation;
            // The conditional branching is short-curcuitting. Thus, self.nodes[dist_idx]
            // is never evalueated when self.nodes.len() <= dist_idx holds.
            let insertion_transition =
                if self.nodes.len() <= dist_idx || self.nodes[dist_idx].has_edge() {
                    prev[node] * config.p_ins
                        + prev[node + 1] * config.p_extend_ins
                        + prev[node + 2] * config.p_del_to_ins
                } else {
                    prev[node..=node + 2].iter().sum::<f64>()
                };
            let insertion_observation = match dist_idx.cmp(&self.nodes.len()) {
                Ordering::Less => self.nodes[dist_idx].insertion(base),
                Ordering::Greater => 0.25,
                Ordering::Equal => 0.25,
            };
            updates[node + 1] = insertion_transition * insertion_observation;
        }
        let d = Self::sum(updates).recip();
        Self::mul(updates, d);
        for (dist_idx, src_nodes) in edges.iter().enumerate() {
            updates[3 * dist_idx + 2] = src_nodes
                .iter()
                .map(|&(src, weight)| {
                    let trans = updates[3 * src] * config.p_del
                        + updates[3 * src + 2] * config.p_extend_del;
                    trans * weight
                })
                .sum::<f64>();
        }
        let c = Self::sum(&updates).recip();
        Self::mul(updates, c);
        (c, d)
    }
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        if self.nodes.is_empty() {
            return DEFAULT_LK;
        }
        // Alignemnts: [mat, ins, del,  mat, ins, del,  ....]
        let mut prev: Vec<f64> = vec![0.; 3 * (self.nodes.len() + 2)];
        // Start from the inserion state.
        let ins_and_del = config.p_ins + config.p_del;
        prev[3 * self.nodes.len() + 4] = config.p_ins / ins_and_del;
        prev[3 * self.nodes.len() + 5] = config.p_del / ins_and_del;
        if prev.len() % 4 != 0 {
            prev.extend(std::iter::repeat(0.).take(4 - prev.len() % 4));
        }
        let mut updated = vec![0.; prev.len()];
        let total = self
            .nodes
            .iter()
            .filter(|n| n.is_head)
            .map(|n| n.weight())
            .sum::<f64>();
        let base_freq = [0.25; 4];
        let edges = {
            let mut edges: Vec<Vec<_>> = vec![vec![]; self.nodes.len() + 2];
            for (from, n) in self.nodes.iter().enumerate() {
                for (&to, &w) in n.edges.iter().zip(n.weights().iter()) {
                    edges[to].push((from, w));
                }
                if n.is_head {
                    edges[from].push((self.nodes.len(), n.weight() / total));
                }
            }
            edges[self.nodes.len()].push((self.nodes.len() + 1, 1.));
            edges
        };
        obs.iter()
            .enumerate()
            .map(|(idx, &base)| {
                updated.iter_mut().for_each(|e| *e = 0.);
                let (c, d) = self.update_row(&mut updated, &prev, base, config, &edges, &base_freq);
                std::mem::swap(&mut prev, &mut updated);
                assert!(c * d > 0.99, "{},{},{},{}", idx, c, d, c * d);
                if idx < obs.len() - 1 {
                    -(c * d).ln()
                } else {
                    -d.ln()
                }
            })
            .sum::<f64>()
    }
}
