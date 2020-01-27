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
    let mut rng: StdRng = SeedableRng::seed_from_u64(121332983);
    let chain_len = 40;
    let template: Vec<_> = (0..40).map(|_| generate_seq(&mut rng, len)).collect();
    let p = Profile {
        sub: 0.0001,
        ins: 0.0001,
        del: 0.0001,
    };
    let template2 = template
        .iter()
        .map(|t| introduce_randomness(t, &mut rng, &p))
        .collect();
    let coverage = 50;
    let mut gen = |ts: &Vec<Vec<u8>>| {
        ts.iter()
            .map(|t| introduce_randomness(&t, &mut rng, &PROFILE))
            .collect::<Vec<_>>()
    };
    let data1: Vec<_> = (0..coverage).map(|_| gen(&template)).collect();
    let data2: Vec<_> = (0..coverage).map(|_| gen(&template2)).collect();
    let data = {
        let mut data = vec![vec![]; chain_len];
        for d in data1.iter().chain(data2.iter()) {
            for (idx, c) in d.iter().enumerate() {
                data[idx].push(c.as_slice());
            }
        }
        data
    };
    let d = template
        .iter()
        .zip(template2.iter())
        .map(|(t1, t2)| edlib_sys::global_dist(t1, t2))
        .enumerate()
        .map(|(idx, d)| {
            println!("{}:{}", idx, d);
            d
        })
        .sum::<u32>();
    println!("{}", d);
    let weight1 = vec![vec![1.; coverage], vec![0.0; coverage]].concat();
    let weight2 = vec![vec![0.1; coverage], vec![0.8; coverage]].concat();
    let mut f = Factory::new();
    let mut buf = vec![];
    let m1: Vec<_> = data
        .iter()
        .map(|chunks| {
            let mut m = f.generate_with_weight_prior(chunks, &weight1, k, &mut buf);
            m.check(&chunks, &DEFAULT_CONFIG, -140.);
            m
        })
        .collect();
    let m2: Vec<_> = data
        .iter()
        .map(|chunks| {
            let mut m = f.generate_with_weight_prior(chunks, &weight2, k, &mut buf);
            m.check(&chunks, &DEFAULT_CONFIG, -140.);
            m
        })
        .collect();
    for (idx, (m1, m2)) in m1.iter().zip(m2.iter()).enumerate() {
        println!("{}\t{}\t{}", idx, m1, m2);
    }
    let num = 500;
    let test: Vec<_> = (0..2 * num)
        .map(|i| {
            if i % 2 == 0 {
                gen(&template)
            } else {
                gen(&template2)
            }
        })
        .collect();
    fn calc_lk(ms: &[DBGHMM], q: &Vec<Vec<u8>>) -> f64 {
        let (lk, len) = ms
            .iter()
            .zip(q.iter())
            .enumerate()
            .filter_map(|(idx, (m, q))| {
                if m.is_broken() {
                    None
                } else {
                    let lk = m.forward(&q, &dbg_hmm::DEFAULT_CONFIG);
                    eprintln!("LK\t{}\t{}", idx, lk);
                    Some(lk)
                }
            })
            .fold((0., 0.), |(lk, num), x| (lk + x, num + 1.));
        lk / len * q.len() as f64
    };
    let correct1 = test
        .par_iter()
        .enumerate()
        .filter(|(idx, _)| idx % 2 == 0)
        .filter(|(_, q)| {
            let f1 = calc_lk(&m1, q);
            let f2 = calc_lk(&m2, q);
            f1 > f2
        })
        .count();
    let correct2 = test
        .par_iter()
        .enumerate()
        .filter(|(idx, _)| idx % 2 == 1)
        .filter(|(_, q)| {
            let f1 = calc_lk(&m1, q);
            let f2 = calc_lk(&m2, q);
            f1 < f2
        })
        .count();
    println!("Fin:{}:{}", correct1, correct2);
}
