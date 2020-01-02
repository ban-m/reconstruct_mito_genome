extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let k = 6;
    let mut rng: StdRng = SeedableRng::seed_from_u64(899_902);
    let p = Profile {
        ins: 0.02,
        del: 0.02,
        sub: 0.02,
    };
    for _ in 0..20 {
        let template1 = generate_seq(&mut rng, len);
        let template2 = gen_sample::introduce_randomness(&template1, &mut rng, &p);
        use std::collections::HashSet;
        let kmer1 = template1
            .windows(k)
            .map(|e| e.to_vec())
            .collect::<HashSet<_>>()
            .len();
        let kmer2 = template2
            .windows(k)
            .map(|e| e.to_vec())
            .collect::<HashSet<_>>()
            .len();
        eprintln!("Dist:{}", edlib_sys::global_dist(&template1, &template2));
        let unit = 20;
        let (cov1, cov2) = (unit, 3 * unit);
        let weight1 = vec![vec![0.8; cov1], vec![0.2; cov2]].concat();
        let weight2 = vec![vec![0.2; cov1], vec![0.8; cov2]].concat();
        let data1: Vec<Vec<_>> = (0..cov1)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let data2: Vec<Vec<_>> = (0..cov2)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let dataset: Vec<_> = data1
            .iter()
            .chain(data2.iter())
            .map(|e| e.as_slice())
            .collect();
        let m1 = DBGHMM::new_with_weight_prior(&dataset, &weight1, k);
        let m2 = DBGHMM::new_with_weight_prior(&dataset, &weight2, k);
        eprintln!("({}){}\t({}){}", kmer1, m1, kmer2, m2);
    }
}
