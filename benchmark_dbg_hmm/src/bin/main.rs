extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let num_seq = 50;
    let mut rng: StdRng = SeedableRng::seed_from_u64(12121899892);
    let template1: Vec<_> = generate_seq(&mut rng, len);
    let p = Profile {
        sub: 0.005,
        ins: 0.005,
        del: 0.005,
    };
    let template2 = introduce_randomness(&template1, &mut rng, &p);
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
        PROFILE.sub, PROFILE.ins, PROFILE.del
    );
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    // for line in &data1 {
    //     println!("Model1\t{}", String::from_utf8_lossy(line));
    // }
    // for line in &data2 {
    //     println!("Model2\t{}", String::from_utf8_lossy(line));
    // }
    let k = 7;
    println!("K={}", k);
    let model1 = DBGHMM::new(&data1, k);
    let model2 = DBGHMM::new(&data2, k);
    let likelihood1 = model1.forward(&template1, &DEFAULT_CONFIG);
    let likelihood2 = model1.forward(&template2, &DEFAULT_CONFIG);
    println!("Model1's likelihood:{:.4}\t{:.4}", likelihood1, likelihood2);
    let likelihood1 = model2.forward(&template1, &DEFAULT_CONFIG);
    let likelihood2 = model2.forward(&template2, &DEFAULT_CONFIG);
    println!("Model2's likelihood:{:.4}\t{:.4}", likelihood1, likelihood2);
    println!("Cross Validation");
    println!("Answer\tModel1\tModel2");
    let cv = cross_validation(&data1, &data2);
    let (mut tot, mut correct) = (0, 0);
    for (ans, m1, m2) in cv {
        let is_ok = (ans == 1 && m1 > m2) || (ans == 2 && m1 < m2);
        tot += 1;
        correct += is_ok as u32;
        let sign = if is_ok { "OK" } else { "NG" };
        println!("{}\t{:.4}\t{:.4}\t{}", ans, m1, m2, sign);
    }
    println!("{} out of {}", correct, tot);
    let cv = cross_validation_nn(&data1, &data2);
    let (mut tot, mut correct) = (0, 0);
    for (ans, m1, m2) in cv {
        let is_ok = (ans == 1 && m1 > m2) || (ans == 2 && m1 < m2);
        tot += 1;
        correct += is_ok as u32;
        let sign = if is_ok { "OK" } else { "NG" };
        println!("{}\t{:.4}\t{:.4}\t{}", ans, m1, m2, sign);
    }
    println!("{} out of {}", correct, tot);
}

// If the query's nearest neighbor is data1, (1,0) else (0,1).
fn cross_validation_nn(data1: &[Vec<u8>], data2: &[Vec<u8>]) -> Vec<(usize, u8, u8)> {
    let len = data1.len();
    let mut res = vec![];
    for i in 0..len {
        let test1 = &data1[i];
        let test2 = &data2[i];
        let (nn_from1_test1, nn_from1_test2) = data1
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx != i)
            .map(|(_, e)| e)
            .map(|seq| {
                (
                    edlib_sys::global_dist(seq, test1),
                    edlib_sys::global_dist(seq, test2),
                )
            })
            .fold((160, 160), |(one, two), (x, y)| (one.min(x), two.min(y)));
        let (nn_from2_test1, nn_from2_test2) = data2
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx != i)
            .map(|(_, e)| e)
            .map(|seq| {
                (
                    edlib_sys::global_dist(seq, test1),
                    edlib_sys::global_dist(seq, test2),
                )
            })
            .fold((160, 160), |(one, two), (x, y)| (one.min(x), two.min(y)));
        if nn_from1_test1 < nn_from2_test1 {
            res.push((1, 1, 0));
        } else {
            res.push((1, 0, 1));
        }
        if nn_from1_test2 < nn_from2_test2 {
            res.push((2, 1, 0));
        } else {
            res.push((2, 0, 1));
        }
    }
    res.sort_by_key(|e| e.0);
    res
}

fn cross_validation(data1: &[Vec<u8>], data2: &[Vec<u8>]) -> Vec<(usize, f64, f64)> {
    let len = data1.len();
    let mut res = vec![];
    let k = 7;
    for i in 0..len {
        let test1 = &data1[i];
        let test2 = &data2[i];
        let data1: Vec<_> = data1
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx != i)
            .map(|(_, e)| e)
            .cloned()
            .collect();
        let data2: Vec<_> = data2
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx != i)
            .map(|(_, e)| e)
            .cloned()
            .collect();
        let m1 = DBGHMM::new(&data1, k);
        let m2 = DBGHMM::new(&data2, k);
        let m1_for_1 = m1.forward(test1, &DEFAULT_CONFIG);
        let m1_for_2 = m1.forward(test2, &DEFAULT_CONFIG);
        let m2_for_1 = m2.forward(test1, &DEFAULT_CONFIG);
        let m2_for_2 = m2.forward(test2, &DEFAULT_CONFIG);
        res.push((1, m1_for_1, m2_for_1));
        res.push((2, m1_for_2, m2_for_2));
    }
    res.sort_by_key(|e| e.0);
    res
}
