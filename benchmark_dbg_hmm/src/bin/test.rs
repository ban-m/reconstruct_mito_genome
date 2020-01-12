extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    let len = 150;
    let k = 6;
    // let mut rng: StdRng = SeedableRng::seed_from_u64(12132334);
    let mut rng: StdRng = SeedableRng::seed_from_u64(121332983);
    let template = generate_seq(&mut rng, len);
    let unit = 40;
    let data1: Vec<_> = (0..unit)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let weight1 = vec![1.; unit];
    let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
    let m1 = DBGHMM::new_with_weight_prior(&data1, &weight1, k);
    m1.dump(&template);
    // eprintln!("{}", m1);
    for d in &data1 {
        eprintln!("{}", String::from_utf8_lossy(d));
    }
    let num = 20;
    let ave = (0..num)
        .map(|_| {
            // let d = introduce_errors(&template, &mut rng, 4, 4, 4);
            let d = introduce_randomness(&template, &mut rng, &PROFILE);
            let f = m1.forward(&d, &dbg_hmm::DEFAULT_CONFIG);
            eprintln!("{:.3}", f);
            f
        })
        .sum::<f64>()
        / num as f64;
    eprintln!("Ave:{}", ave);
}
