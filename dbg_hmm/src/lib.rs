#![feature(test)]
//! A tiny implementation for de Bruijn graph with Hidden Markov model.
//! Currently, this implementation is minimal. In other words, it exposes only one struct with just two methods:
//! [DeBruijnGraphHiddenMarkovModel] -- Yes, it is too long -- ,[constructor](DeBruijnGraphHiddenMarkovModel::new),
//! and the [forward algorithm](DeBruijnGraphHiddenMarkovModel::forward)
//! to calculate the probability this graph would generate the given observation.
//! As a shorthand for the vary long name, I also supply [DBGHMM] as a alias for [DeBruijnGraphHiddenMarkovModel].
extern crate edlib_sys;
extern crate env_logger;
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
mod kmer;
use kmer::Kmer;
// This setting is determined by experimentally.
mod config;
pub use config::*;
use packed_simd::f64x4 as f64s;
pub type DBGHMM = DeBruijnGraphHiddenMarkovModel;
mod factory;
pub use factory::Factory;

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
    fn is_non_zero(idx: usize, xs: &[f64]) -> bool {
        xs[idx * 4..(idx + 1) * 4].iter().any(|&e| e > 0.00001)
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
    fn initialize(&self, tip: &[u8], config: &Config) -> (f64, f64, Vec<f64>) {
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
        let hat: Vec<_> = tilde.flat_map(|init| vec![init, 0., 0., 0.]).collect();
        let minus_ln_c = 0.;
        (minus_ln_c, minus_ln_d, hat)
    }
    /// Return the total weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    #[cfg(target_feature = "sse")]
    pub fn forward(&self, obs: &[u8], config: &Config) -> f64 {
        assert!(obs.len() > self.k);
        if self.weight < WEIGHT_THR {
            return LOW_LIKELIHOOD;
        }
        // Alignemnts: [mat, ins, del, 0, mat, ins, del, 0, ....]
        let (mut cs, mut ds, mut prev) = self.initialize(&obs[..self.k], config);
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::Xoroshiro128StarStar;
    #[test]
    fn works() {}
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
