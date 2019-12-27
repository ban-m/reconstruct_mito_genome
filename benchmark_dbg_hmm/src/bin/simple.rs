extern crate dbg_hmm;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let num_seq = 20;
    let k = 6;
    let mut rng: StdRng = SeedableRng::seed_from_u64(899_892_321);
    let template1: Vec<_> = generate_seq(&mut rng, len);
    eprintln!("{}", String::from_utf8_lossy(&template1));
    let template2 = gen_sample::introduce_errors(&template1, &mut rng, 0, 0, 1);
    let data1: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let m1 = {
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let weight = vec![1.; num_seq];
        DBGHMM::new_with_weight_prior(&data1, &weight, k)
    };
    let data2: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    let m2 = {
        let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
        let weight = vec![1.; num_seq];
        DBGHMM::new_with_weight_prior(&data2, &weight, k)
    };
    eprintln!("{}\t{}", m1, m2);
    assert!(m1.is_connected(), "M1");
    assert!(m2.is_connected(), "M2");
    let correct = (0..100)
        .filter(|e| {
            eprintln!("------");
            if e % 2 == 0 {
                let q = gen_sample::introduce_randomness(&template1, &mut rng, &PROFILE);
                let s1 = m1.forward(&q, &DEFAULT_CONFIG);
                let s2 = m2.forward(&q, &DEFAULT_CONFIG);
                eprintln!("{}\t{}\t{:6.2}\t{:6.2}", e % 2, s1 > s2, s1, s2);
                viterbi_print(&q, &m1);
                viterbi_print(&q, &m2);
                s1 > s2
            } else {
                let q = gen_sample::introduce_randomness(&template2, &mut rng, &PROFILE);
                let s1 = m1.forward(&q, &DEFAULT_CONFIG);
                let s2 = m2.forward(&q, &DEFAULT_CONFIG);
                eprintln!("{}\t{}\t{:6.2}\t{:6.2}", e % 2, s1 > s2, s1, s2);
                viterbi_print(&q, &m1);
                viterbi_print(&q, &m2);
                s1 < s2
            }
        })
        .count();
    eprintln!("{}", correct as f64 / 100.);
}

fn viterbi_print(q: &[u8], m: &DBGHMM) {
    let s = m.forward_exp(&q, &DEFAULT_CONFIG);
    let (max, viterbi) = m.viterbi(&q, &DEFAULT_CONFIG);
    eprintln!("Forward:{} vs Viterbi{}", s, max);
    let mut pos = 0;
    let (mut refr, mut query) = (String::new(), String::new());
    for (base, state) in viterbi {
        match state {
            0 => {
                refr.push(base as char);
                query.push(q[pos] as char);
                pos += 1;
            }
            1 => {
                refr.push('-');
                query.push(q[pos] as char);
                pos += 1;
            }
            2 => {
                refr.push(base as char);
                query.push('-');
            }
            _ => {}
        }
    }
    eprintln!("{}", query);
    eprintln!("{}", refr);
}
