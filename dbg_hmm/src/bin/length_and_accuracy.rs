extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate log;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
fn main() {
    //env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let min_len = 50;
    let max_len = 300;
    let by = 5;
    let num_seq = 10;
    let test_num = 1000;
    let k = 6;
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12218993492);
    println!("Length\tDist\tProposed\tNaive\tErrors");
    let p = gen_sample::Profile {
        sub: 0.003,
        ins: 0.003,
        del: 0.003,
    };
    (0..100).for_each(|_| {
        (min_len / by..max_len / by)
            .map(|e| simulate(e * by, num_seq, test_num, k, &mut rng, &p))
            .for_each(|(len, dist, acc1, acc2)| {
                println!("{}\t{}\t{}\t{}\t0.9", len, dist, acc1, acc2)
            })
    });
    let p = gen_sample::Profile {
        sub: 0.01,
        ins: 0.01,
        del: 0.01,
    };
    (0..100).for_each(|_| {
        (min_len / by..max_len / by)
            .map(|e| simulate(e * by, num_seq, test_num, k, &mut rng, &p))
            .for_each(|(len, dist, acc1, acc2)| {
                println!("{}\t{}\t{}\t{}\t0.9", len, dist, acc1, acc2)
            })
    });
}

fn simulate<T: Rng>(
    len: usize,
    num_seq: usize,
    test_num: usize,
    k: usize,
    rng: &mut T,
    p: &gen_sample::Profile,
) -> (usize, usize, f64, f64) {
    let template1 = dbg_hmm::gen_sample::generate_seq(rng, len);
    let template2 = dbg_hmm::gen_sample::introduce_randomness(&template1, rng, p);
    let dist = edlib_sys::global(&template1, &template2)
        .into_iter()
        .filter(|&e| e != 0)
        .count();
    let p = &gen_sample::PROFILE;
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, rng, p))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, rng, p))
        .collect();
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
            let l1 = model1.forward(&test, &DEFAULT_CONFIG);
            let l2 = model2.forward(&test, &DEFAULT_CONFIG);
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
