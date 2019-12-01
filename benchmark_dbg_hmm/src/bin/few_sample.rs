extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
use std::time;
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let len = 150;
    let num_seq = (5..15).collect::<Vec<usize>>();
    let test_num = 1000;
    let k = 6;
    let sample_num: Vec<u64> = (0..800).collect();
    let p = &gen_sample::Profile {
        sub: 0.003,
        ins: 0.004,
        del: 0.003,
    };
    let c = if args[1] == "default" {
        DEFAULT_CONFIG
    } else {
        PACBIO_CONFIG
    };
    eprintln!("Introducing errors");
    eprintln!("Sub:{}\tIns:{}\tDel:{}", p.sub, p.ins, p.del);
    let s = &gen_sample::PROFILE;
    eprintln!("Sequencing errors");
    eprintln!("Sub:{}\tIns:{}\tDel:{}", s.sub, s.ins, s.del);
    eprintln!("K={}", k);
    let params: Vec<_> = sample_num
        .iter()
        .flat_map(|&e| num_seq.iter().map(|&n| (e, n)).collect::<Vec<_>>())
        .collect();
    let result: Vec<_> = params
        .par_iter()
        .map(|&(seed, num_seq)| benchmark(p, s, seed, k, len, num_seq, test_num, &c))
        .collect();
    println!("HMM\tAln\tDist\tCoverage");
    for (hmm, aln, dist, num_seq) in result {
        println!("{}\t{}\t{}\t{}", hmm, aln, dist, num_seq);
    }
}
fn benchmark(
    p: &gen_sample::Profile,
    s: &gen_sample::Profile,
    seed: u64,
    k: usize,
    len: usize,
    num_seq: usize,
    test_num: usize,
    config: &Config,
) -> (f64, f64, u32, usize) {
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(122_492 + seed);
    let template1 = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let template2 = dbg_hmm::gen_sample::introduce_randomness(&template1, &mut rng, &p);
    let dist = edlib_sys::global_dist(&template1, &template2);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, &mut rng, s))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, &mut rng, s))
        .collect();
    let model1 = DBGHMM::new(&data1, k);
    let model2 = DBGHMM::new(&data2, k);
    let tests: Vec<_> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(0.5) {
                (1, gen_sample::introduce_randomness(&template1, &mut rng, p))
            } else {
                (2, gen_sample::introduce_randomness(&template2, &mut rng, p))
            }
        })
        .collect();
    let start = time::Instant::now();
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = model1.forward(&test, config);
            let l2 = model2.forward(&test, config);
            (l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2)
        })
        .count();
    let time = time::Instant::now() - start;
    let time_par_case = time.as_millis() as f64 / test_num as f64;
    eprintln!("{:?}({:.3}millis/sample)", time, time_par_case);
    let hmm = correct as f64 / test_num as f64;
    let correct = tests
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
            let is_1 = match from1.cmp(&from2) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Equal => rng.gen_bool(0.5),
                std::cmp::Ordering::Greater => false,
            };
            (is_1 && *ans == 1) || (!is_1 && *ans == 2)
        })
        .count();
    let aln = correct as f64 / test_num as f64;
    (hmm, aln, dist, num_seq)
}
