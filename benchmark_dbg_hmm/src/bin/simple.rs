extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let len = 150;
    let k = 6;
    let chain = 40;
    let seed: u64 = std::env::args()
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap_or(899_902);
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
    let p = Profile {
        ins: 0.0005,
        del: 0.0005,
        sub: 0.0005,
    };
    let template1: Vec<_> = (0..chain).map(|_| generate_seq(&mut rng, len)).collect();
    let template2: Vec<_> = template1
        .iter()
        .map(|seq| gen_sample::introduce_randomness(seq, &mut rng, &p))
        .collect();
    use std::collections::HashSet;
    let kmer1: Vec<_> = template1
        .iter()
        .map(|e| {
            e.windows(k)
                .map(|e| e.to_vec())
                .collect::<HashSet<_>>()
                .len()
        })
        .collect();
    let kmer2: Vec<_> = template2
        .iter()
        .map(|e| {
            e.windows(k)
                .map(|e| e.to_vec())
                .collect::<HashSet<_>>()
                .len()
        })
        .collect();
    let dist = template1
        .iter()
        .zip(template2.iter())
        .map(|(s, q)| edlib_sys::global_dist(s, q))
        .sum::<u32>();
    eprintln!("Dist:{}", dist);
    let unit = 30;
    let (cov1, cov2) = (unit, 10 * unit);
    let data1: Vec<Vec<_>> = (0..cov1)
        .map(|_| {
            template1
                .iter()
                .map(|seq| introduce_randomness(seq, &mut rng, &PROFILE))
                .collect()
        })
        .collect();
    let data2: Vec<Vec<_>> = (0..cov2)
        .map(|_| {
            template2
                .iter()
                .map(|seq| introduce_randomness(seq, &mut rng, &PROFILE))
                .collect()
        })
        .collect();
    let weight1 = vec![1.; cov1];
    let m1: Vec<_> = (0..chain)
        .map(|c| {
            let data: Vec<_> = data1.iter().map(|e| e[c].as_slice()).collect();
            DBGHMM::new_with_weight_prior(&data, &weight1, k)
        })
        .collect();
    let weight2 = vec![1.; cov2];
    let m2: Vec<_> = (0..chain)
        .map(|c| {
            let data: Vec<_> = data2.iter().map(|e| e[c].as_slice()).collect();
            DBGHMM::new_with_weight_prior(&data, &weight2, k)
        })
        .collect();
    for ((m1, k1), (m2, k2)) in m1.iter().zip(kmer1.iter()).zip(m2.iter().zip(kmer2.iter())) {
        eprintln!("({}){}\t({}){}", k1, m1, k2, m2);
    }
    fn offset(x: f64) -> f64 {
        const A_PRIOR: f64 = -0.1;
        const B_PRIOR: f64 = 3.1;
        // const A_PRIOR: f64 = -0.089;
        // const B_PRIOR: f64 = 3.1;
        // const A_PRIOR: f64 = -0.0165;
        // const B_PRIOR: f64 = 0.7;
        // const A_PRIOR: f64 = -0.003;
        // const B_PRIOR: f64 = 1.10;
        // const A_PRIOR: f64 = -0.0089;
        // const B_PRIOR: f64 = 2.08;
        (B_PRIOR + x * A_PRIOR).exp()
    }
    eprintln!("{}\t{}", cov1, cov2);
    let sum1 = data1
        .par_iter()
        .map(|read| {
            let lk1 = read
                .iter()
                .zip(m1.iter())
                .map(|(s, m)| m.forward(s, &DEFAULT_CONFIG) + offset(m.weight()))
                .sum::<f64>();
            let lk2 = read
                .iter()
                .zip(m2.iter())
                .map(|(s, m)| m.forward(s, &DEFAULT_CONFIG) + offset(m.weight()))
                .sum::<f64>();
            eprintln!("{}\t{}\t{}", lk1, lk2, lk1 > lk2);
            (lk1 > lk2) as usize
        })
        .sum::<usize>();
    eprintln!("========");
    let sum2 = data2
        .par_iter()
        .map(|read| {
            let lk1 = read
                .iter()
                .zip(m1.iter())
                .map(|(s, m)| m.forward(s, &DEFAULT_CONFIG) + offset(m.weight()))
                .sum::<f64>();
            let lk2 = read
                .iter()
                .zip(m2.iter())
                .map(|(s, m)| m.forward(s, &DEFAULT_CONFIG) + offset(m.weight()))
                .sum::<f64>();
            eprintln!("{}\t{}\t{}", lk1, lk2, lk1 < lk2);
            (lk1 < lk2) as usize
        })
        .sum::<usize>();
    println!("{}:{}:{}/{},{}/{}", seed, dist, sum1, cov1, sum2, cov2);
}
