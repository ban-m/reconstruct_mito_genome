extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::gen_sample::*;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
// fn alignment<F>(xs: &[u8], ys: &[u8], del: f64, ins: f64, score: F) -> f64
// where
//     F: Fn(u8, u8) -> f64,
// {
//     let mut dp = vec![vec![0.; ys.len() + 1]; xs.len() + 1];
//     for i in 0..xs.len() {
//         dp[i + 1][0] = ins * (i + 1) as f64;
//     }
//     // for j in 0..ys.len() {
//     //     dp[0][j + 1] = del * (j + 1) as f64;
//     // }
//     for i in 0..xs.len() {
//         for j in 0..ys.len() {
//             dp[i + 1][j + 1] = (dp[i][j] + score(xs[i], ys[j]))
//                 .max(dp[i + 1][j] + del)
//                 .max(dp[i][j + 1] + ins);
//         }
//     }
//     *dp[xs.len()]
//         .iter()
//         .max_by(|a, b| a.partial_cmp(&b).unwrap())
//         .unwrap()
// }

fn main() {
    let s = &PROFILE;
    let p = &Profile {
        sub: 0.003,
        ins: 0.003,
        del: 0.003,
    };
    let num_seq = 50;
    let len = 150;
    let test_num = 1000;
    let config = &DEFAULT_CONFIG;
    let seeds: Vec<_> = (0..1000).collect();
    //seeds.par_iter().for_each(|&i| {
    let i = 460;
    //    let i = 10;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(i);
    let template1 = gen_sample::generate_seq(&mut rng, len);
    let template2 = gen_sample::introduce_randomness(&template1, &mut rng, &p);
    let dist = edlib_sys::global_dist(&template1, &template2);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, &mut rng, s))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, &mut rng, s))
        .collect();
    // let data1: Vec<Vec<_>> = data1
    //     .into_iter()
    //     .map(|e| e.into_iter().take(50).collect())
    //     .collect();
    // let data2: Vec<Vec<_>> = data2
    //     .into_iter()
    //     .map(|e| e.into_iter().take(50).collect())
    //     .collect();
    let model1 = POA::generate_vec(&data1);
    let model2 = POA::generate_vec(&data2);
    let tests: Vec<_> = (0..test_num)
        .map(|i| {
            let is_one = i < test_num / 2;
            let q = if is_one {
                gen_sample::introduce_randomness(&template1, &mut rng, p)
            } else {
                gen_sample::introduce_randomness(&template2, &mut rng, p)
            };
            // let q: Vec<_> = q.into_iter().take(50).collect();
            (if is_one { 1 } else { 2 }, q)
        })
        .collect();
    let mut as_one = 0;
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = model1.forward(&test, config);
            let l2 = model2.forward(&test, config);
            if l2 < l1 {
                as_one += 1;
            }
            // eprintln!("{}\t{:.3}\t{:.3}", ans, l1, l2);
            (l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2)
        })
        .count();
    let hmm = correct as f64 / test_num as f64;
    if hmm < 0.6 && dist == 1 {
        eprintln!("Acc:{}\t{}\t{}\n{}\n{}", as_one, hmm, i, model1, model2);
        eprintln!("REF:{}", String::from_utf8_lossy(&template1));
        eprintln!("REF:{}", String::from_utf8_lossy(&template2));
        eprintln!(
            "MOD:{}",
            String::from_utf8_lossy(&model1.consensus().unwrap())
        );
        eprintln!(
            "MOD:{}",
            String::from_utf8_lossy(&model2.consensus().unwrap())
        );
    }
    //    });
    eprintln!("{:?}", model1);
    // for d in data2 {
    //     eprintln!("{}", String::from_utf8_lossy(&d));
    // }
}
