extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use std::time;
fn main() {
    let len = 100;
    let num_seq = 10;
    let test_num = 1000;
    let k = 6;
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12218993492);
    let template1 = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let p = gen_sample::Profile {
        sub: 0.003,
        ins: 0.004,
        del: 0.003,
    };
    let template2 = dbg_hmm::gen_sample::introduce_randomness(&template1, &mut rng, &p);
    println!("Introducing errors");
    println!("Sub:{}\tIns:{}\tDel:{}", p.sub, p.ins, p.del);
    println!("Template1\t{}", String::from_utf8_lossy(&template1));
    println!("Template2\t{}", String::from_utf8_lossy(&template2));
    println!(
        "Distance:{}",
        edlib_sys::global(&template1, &template2)
            .into_iter()
            .filter(|&e| e != 0)
            .count()
    );
    // let p = &dbg_hmm::gen_sample::Profile {
    //     sub: 0.03,
    //     del: 0.05,
    //     ins: 0.06,
    // };
    let p = &gen_sample::PROFILE;
    println!("Sequencing errors");
    println!("Sub:{}\tIns:{}\tDel:{}", p.sub, p.ins, p.del);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, &mut rng, p))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, &mut rng, p))
        .collect();
    println!("K={}", k);
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
            let l1 = model1.forward(&test, &DEFAULT_CONFIG);
            let l2 = model2.forward(&test, &DEFAULT_CONFIG);
            // eprintln!("{}\t{:.1}\t{:.1}", ans, l1, l2);
            (l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2)
        })
        .count();
    let time = time::Instant::now() - start;
    println!("Answer\tModel1\tModel2");
    println!(
        "{} out of {} ({:.3} %). {:?}({:.3}millis/sample)",
        correct,
        test_num,
        correct as f64 / test_num as f64,
        time,
        time.as_millis() as f64 / test_num as f64,
    );
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
            let is_1 = if from1 < from2 {
                true
            }else if from1 == from2 {
                rng.gen_bool(0.5)
            }else{
                false
            };
            (is_1 && *ans == 1) || (!is_1 && *ans == 2)
        })
        .count();
    println!(
        "{} out of {} ({} %)",
        correct,
        test_num,
        correct as f64 / test_num as f64
    );
}
