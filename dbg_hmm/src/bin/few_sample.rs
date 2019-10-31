extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate log;
extern crate rand;
use dbg_hmm::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
fn main() {
    //env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 25;
    let test_num = 1000;
    let k = 6;
    let mut rng: StdRng = SeedableRng::seed_from_u64(12121899892);
    let template1 = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
    let p = gen_sample::Profile {
        sub: 0.005,
        ins: 0.005,
        del: 0.005,
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
    println!("Sequencing errors");
    println!(
        "Sub:{}\tIns:{}\tDel:{}",
        gen_sample::PROFILE.sub,
        gen_sample::PROFILE.ins,
        gen_sample::PROFILE.del
    );
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| gen_sample::introduce_randomness(&template2, &mut rng, &gen_sample::PROFILE))
        .collect();
   println!("K={}", k);
    let model1 = DBGHMM::new(&data1, k);
    let model2 = DBGHMM::new(&data2, k);
    let tests: Vec<_> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(0.5) {
                (
                    1,
                    gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE),
                )
            } else {
                (
                    2,
                    gen_sample::introduce_randomness(&template2, &mut rng, &gen_sample::PROFILE),
                )
            }
        })
        .collect();
    println!("Answer\tModel1\tModel2");
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = model1.forward(&test, &DEFAULT_CONFIG);
            let l2 = model2.forward(&test, &DEFAULT_CONFIG);
            (l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2)
        })
        .count();
    println!(
        "{} out of {} ({} %)",
        correct,
        test_num,
        correct as f64 / test_num as f64
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
            (from1 < from2 && *ans == 1) || (from2 < from1 && *ans == 2)
        })
        .count();
    println!(
        "{} out of {} ({} %)",
        correct,
        test_num,
        correct as f64 / test_num as f64
    );
}
