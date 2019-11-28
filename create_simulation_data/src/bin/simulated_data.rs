extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use create_simulation_data::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let k = 6;
    let len = 100;
    let min_coverage = 10;
    let max_coverage = 31;
    let test_num = 100;
    let sample_num: Vec<(u64, usize)> = (0..500)
        .flat_map(|e| {
            (min_coverage..max_coverage)
                .map(|c| (e, c))
                .collect::<Vec<_>>()
        })
        .collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let result: Vec<_> = sample_num
        .into_par_iter()
        .map(|(seed, coverage)| {
            let (hmm, aln, dist) = benchmark(seed, p, coverage, test_num, chain_len, k, len);
            (hmm, aln, dist, coverage)
        })
        .collect();
    println!("HMM\tAln\tDist\tCoverage\tLength");
    for (hmm, aln, dist, coverage) in result {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            hmm,
            aln,
            dist,
            coverage,
            len * chain_len
        );
    }
}
fn benchmark(
    seed: u64,
    p: &gen_sample::Profile,
    coverage: usize,
    test_num: usize,
    chain_len: usize,
    k: usize,
    len: usize,
) -> (f64, f64, u32) {
    let seed = 100342374 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let template1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let template2: Vec<_> = template1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let dist = template1
        .iter()
        .zip(template2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .sum::<u32>();
    let prob_0 = 0.5;
    let (dataset, label, answer, border) =
        generate_dataset(&template1, &template2, coverage, test_num, &mut rng, prob_0);
    let em_pred = em_solve(&dataset, &label, border, k, 10);
    let correct = em_pred
        .iter()
        .zip(answer.iter())
        .filter(|(p, a)| a == p)
        .count();
    let hmm = correct as f64 / test_num as f64;
    let aln = align_pred(&dataset, &label, border);
    let correct = aln
        .iter()
        .zip(answer.iter())
        .filter(|(p, a)| p == a)
        .count();
    let aln = correct as f64 / test_num as f64;
    (hmm, aln, dist)
}

fn align_pred(data: &[Vec<Vec<u8>>], label: &[u8], border: usize) -> Vec<u8> {
    let data = data
        .iter()
        .map(|read| {
            read.iter()
                .flat_map(|e| e.iter().map(|&e| e))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let d0: Vec<&Vec<_>> = data
        .iter()
        .zip(label.iter())
        .filter_map(|(r, &b)| if b == 0 { Some(r) } else { None })
        .collect();
    let d1: Vec<&Vec<_>> = data
        .iter()
        .zip(label.iter())
        .filter_map(|(r, &b)| if b != 0 { Some(r) } else { None })
        .collect();
    let mut assign: Vec<_> = label
        .iter()
        .copied()
        .chain(data.iter().skip(border).map(|query| {
            let min0 = d0.iter().map(|r| edlib_sys::global_dist(r, query)).min();
            let min1 = d1.iter().map(|r| edlib_sys::global_dist(r, query)).min();
            if min0 < min1 {
                0
            } else {
                1
            }
        }))
        .collect();
    let mut updated = true;
    while updated {
        let d0: Vec<&Vec<_>> = data
            .iter()
            .zip(assign.iter())
            .filter_map(|(r, &b)| if b == 0 { Some(r) } else { None })
            .collect();
        let d1: Vec<&Vec<_>> = data
            .iter()
            .zip(assign.iter())
            .filter_map(|(r, &b)| if b != 0 { Some(r) } else { None })
            .collect();
        updated = assign
            .par_iter_mut()
            .enumerate()
            .skip(border)
            .map(|(idx, a)| {
                let query = &data[idx];
                let min0 = d0.iter().map(|r| edlib_sys::global_dist(r, query)).min();
                let min1 = d1.iter().map(|r| edlib_sys::global_dist(r, query)).min();
                let p = if min0 < min1 { 0 } else { 1 };
                let updated = p != *a;
                *a = p;
                updated
            })
            .reduce(|| false, |p, q| q | p);
    }
    assign.into_iter().skip(border).collect()
}
