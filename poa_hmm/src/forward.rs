use crate::Config;
use crate::PartialOrderAlignment;
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
    fn update(
        &self,
        updates: &mut [f64],
        prev: &[f64],
        base: u8,
        config: &Config,
        edges: &[Vec<(usize, f64)>],
    ) -> (f64, f64) {
        debug_assert!((1. - prev.iter().sum::<f64>()).abs() < SMALL);
        for (dist_idx, (dist, froms)) in self.nodes.iter().zip(edges.iter()).enumerate() {
            let node = 3 * dist_idx;
            let del_to_match = 1. - config.p_extend_del - config.p_del_to_ins;
            let (match_state, insertion_state) = froms
                .iter()
                .map(|&(src, weight)| {
                    let src_node = &self.nodes[src];
                    let src = src * 3;
                    let f_dist = prev[src] * config.p_match
                        + prev[src + 1] * (1. - config.p_extend_ins)
                        + prev[src + 2] * del_to_match;
                    // let m = f_dist * dist.prob(base, config) * weight;
                    let m = f_dist * dist.prob_with(base, config, src_node) * weight;
                    let i = prev[node + 2] * config.p_del_to_ins * weight;
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
    #[cfg(target_feature = "sse")]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        // eprintln!("{}", String::from_utf8_lossy(obs));
        // Alignemnts: [mat, ins, del,  mat, ins, del,  ....]
        let mut prev: Vec<f64> = self
            .nodes
            .iter()
            .map(|n| {
                if n.is_head {
                    n.prob(obs[0], config)
                } else {
                    0.
                }
            })
            .flat_map(|x| vec![x, 0., 0.])
            .collect();
        if prev.len() % 4 != 0 {
            prev.extend(std::iter::repeat(0.).take(4 - prev.len() % 4));
        }
        let c = prev.iter().sum::<f64>();
        prev.iter_mut().for_each(|e| *e /= c);
        let mut updated = vec![0.; prev.len()];
        let edges = {
            let mut edges: Vec<Vec<_>> = vec![vec![]; self.nodes.len()];
            for (from, n) in self.nodes.iter().enumerate() {
                for (&to, &w) in n.edges.iter().zip(n.weights().iter()) {
                    edges[to].push((from, w));
                }
            }
            edges
        };
        let lk = obs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(idx, &base)| {
                updated.iter_mut().for_each(|e| *e = 0.);
                let (c, d) = self.update(&mut updated, &prev, base, config, &edges);
                std::mem::swap(&mut prev, &mut updated);
                assert!(c * d > 1.);
                if idx < obs.len() - 1 {
                    -(c * d).ln()
                } else {
                    -d.ln()
                }
            })
            .sum::<f64>();
        c.ln() + lk
    }
}
