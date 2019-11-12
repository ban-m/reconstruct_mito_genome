extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args:Vec<_> = std::env::args().collect();
    let min_len = 50;
    let max_len = 300;
    let by = 5;
    let num_seq = 15;
    let test_num = 500;
    let k = 6;
    let rep = (0..50).collect::<Vec<u64>>();
    let c = if args[1] == "default" {
        DEFAULT_CONFIG
    } else {
        PACBIO_CONFIG
    };
    let dists = vec![
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 1),
        (1, 2, 1),
        (2, 1, 1),
        (1, 1, 2),
        (2, 2, 1),
        (2, 1, 2),
        (1, 2, 2),
        (2, 2, 2),
        (3, 3, 3),
        (4, 4, 4),
    ];
    println!("Length\tDist\tProposed\tNaive");
    let result: Vec<_> = rep
        .par_iter()
        .flat_map(|e| {
            let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12218993492 + e);
            let mut res = vec![];
            for len in min_len / by..max_len / by {
                for &d in &dists {
                    res.push(simulate(len * by, num_seq, test_num, k, &mut rng, d, &c));
                }
            }
            res
        })
        .collect();
    for (len, dist, acc1, acc2) in result {
        println!("{}\t{}\t{}\t{}", len, dist, acc1, acc2)
    }
}

fn simulate<T: Rng>(
    len: usize,
    num_seq: usize,
    test_num: usize,
    k: usize,
    rng: &mut T,
    (sub, del, ins): (usize, usize, usize),
    config: &dbg_hmm::Config,
) -> (usize, usize, f64, f64) {
    let template1 = dbg_hmm::gen_sample::generate_seq(rng, len);
    let template2 = dbg_hmm::gen_sample::introduce_errors(&template1, rng, sub, del, ins);
    let p = &gen_sample::PROFILE;
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, rng, p))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, rng, p))
        .collect();
    let dist = sub + del + ins;
    let model1 = DBGHMM::new(&data1, k);
    let model2 = DBGHMM::new(&data2, k);
    let tests: Vec<_> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(0.5) {
                (1, gen_sample::introduce_randomness(&template1, rng, p))
            } else {
                (2, gen_sample::introduce_randomness(&template2, rng, p))
            }
        })
        .collect();
    let proposed = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = model1.forward(&test, config);
            let l2 = model2.forward(&test, config);
            // eprintln!("{}\t{:.1}\t{:.1}", ans, l1, l2);
            (l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2)
        })
        .count() as f64
        / test_num as f64;
    let naive = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let from1 = data1
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, test))
                .min()
                .unwrap();
            let from2 = data2
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, test))
                .min()
                .unwrap();
            let is_1 = if from1 < from2 {
                true
            } else if from1 == from2 {
                rng.gen_bool(0.5)
            } else {
                false
            };
            (is_1 && *ans == 1) || (!is_1 && *ans == 2)
        })
        .count() as f64
        / test_num as f64;
    (len, dist, proposed, naive)
}
