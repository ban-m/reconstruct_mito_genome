extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
#[macro_use]
extern crate log;
extern crate env_logger;
use create_simulation_data::*;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
fn main() {
    use env_logger::Env;
    env_logger::Builder::from_env(Env::default().default_filter_or("debug"))
        .init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 10;
    let k = 6;
    let len = 150;
    let coverage: Vec<_> = (6..20).step_by(2).collect();
    let test_num: Vec<_> = (100..300).step_by(20).collect();
    let skew: Vec<_> = (1..10).map(|e| e as f64 / 10.).collect();
    let sample_num: Vec<_> = (0..12).collect();
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    println!("TP_EM\tTN_EM\tTP_Naive\tTN_Naive\tCoverage\tTestPos\tTestNeg\tLength");
    let tot = coverage.len() * test_num.len();
    let mut proc = 1;
    let param: Vec<_> = sample_num
        .iter()
        .flat_map(|&e| skew.iter().copied().map(|f| (e, f)).collect::<Vec<_>>())
        .collect();
    for &cov in &coverage {
        for &test in &test_num {
            let res: Vec<_> = param
                .iter()
                .filter_map(|&(seed, sk)| {
                    let (em_tp, em_tn, naive_tp, naive_tn, pos, neg) =
                        benchmark(seed, p, cov, test, chain_len, k, len, sk);
                    if pos == 0 || neg == 0 {
                        None
                    } else {
                        let em_tp = em_tp as f64 / pos as f64;
                        let em_tn = em_tn as f64 / neg as f64;
                        let naive_tp = naive_tp as f64 / pos as f64;
                        let naive_tn = naive_tn as f64 / neg as f64;
                        Some((neg, pos, em_tp, em_tn, naive_tp, naive_tn, sk))
                    }
                })
                .collect();
            for (neg, pos, em_tp, em_tn, naive_tp, naive_tn, sk) in res {
                println!(
                    "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}",
                    em_tp,
                    em_tn,
                    naive_tp,
                    naive_tn,
                    cov,
                    pos,
                    neg,
                    len * chain_len,
                    sk,
                );
            }
            debug!("{}/{}", proc, tot);
            proc += 1;
        }
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
    skew: f64,
) -> (usize, usize, usize, usize, usize, usize) {
    let seed = 100342374 + seed;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let templates1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates2: Vec<_> = templates1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let (dataset, assignment, answer, border) =
        generate_dataset(&templates1, &templates2, coverage, test_num, &mut rng, skew);
    let naive_pred = align_solve(&dataset, &assignment, border);
    let em_pred = em_solve_with(&dataset, &assignment, border, k);
    assert_eq!(em_pred.len(), answer.len());
    let (em_tp, em_tn) = calc_metric(&em_pred, &answer);
    let (naive_tp, naive_tn) = calc_metric(&naive_pred, &answer);
    let (pos, neg) = answer.iter().fold(
        (0, 0),
        |(p, n), &x| if x == 0 { (p + 1, n) } else { (p, n + 1) },
    );
    (em_tp, em_tn, naive_tp, naive_tn, pos, neg)
}

fn calc_metric(pred: &[u8], answer: &[u8]) -> (usize, usize) {
    pred.iter()
        .zip(answer.iter())
        .fold((0, 0), |(tp, tn), (&p, &a)| {
            if p == a && a == 0 {
                (tp + 1, tn)
            } else if p == a {
                (tp, tn + 1)
            } else {
                (tp, tn)
            }
        })
}
