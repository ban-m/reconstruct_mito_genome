extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
fn main() {
    let seed = 100342374;
    let chain_len = 20;
    let k = 6;
    let len = 150;
    let max_num = 30;
    let min_num = 25;
    let test_num = 100;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let p = gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let templates2: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates1: Vec<_> = templates2
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    let prob_1 = 0.5;
    println!("Chain Length:{}", chain_len);
    println!("Unit Length:{}", len);
    println!("Number of test cases:{}", test_num);
    println!("k:{}", k);
    println!("Divergence\tSub:{}\tIns:{}\tDel:{}", p.sub, p.ins, p.del);
    println!(
        "# of training data for template1 ~ Unif({},{})",
        min_num, max_num
    );
    println!(
        "# of training data for template2 ~ Unif({},{})",
        min_num, max_num
    );
    println!("Template1:Template2={}:{}",prob_1,1.-prob_1);
    println!("SegmentID\tDivergence");
    for (idx, (t1, t2)) in templates1.iter().zip(templates2.iter()).enumerate() {
        println!("{}\t{}", idx, edlib_sys::global_dist(t1, t2));
        eprintln!("{}\t{}", idx, edlib_sys::global_dist(t1, t2));
    }
    let gen_sample::Profile { sub, ins, del } = gen_sample::PROFILE;
    println!("SeqErrors\tSub:{}\tIns:{}\tDel:{}", sub, ins, del);
    let (dataset2, model2) = generate_dataset(&templates2, min_num, max_num, &mut rng, k);
    let (dataset1, model1) = generate_dataset(&templates1, min_num, max_num, &mut rng, k);
    let tests: Vec<(_, Vec<_>)> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(prob_1) {
                let d = templates1
                    .iter()
                    .map(|e| gen_sample::introduce_randomness(e, &mut rng, &gen_sample::PROFILE))
                    .collect();
                (1, d)
            } else {
                let d = templates2
                    .iter()
                    .map(|e| gen_sample::introduce_randomness(e, &mut rng, &gen_sample::PROFILE))
                    .collect();
                (2, d)
            }
        })
        .collect();
    let start = std::time::Instant::now();
    println!("answer\tpredict\tw1\tw2");
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = predict(&model1, test);
            let l2 = predict(&model2, test);
            let (l1, l2) = merge_predict(&l1, &l2);
            let p = if l2 < l1 { 1 } else { 2 };
            println!("{}\t{}\t{:.4}\t{:.4}", ans, p, l1, l2);
            *ans == p
        })
        .count();
    let ela = std::time::Instant::now() - start;
    println!(
        "Elapsed Time:{:?}/{}\t{:.2}millis/case",
        ela,
        tests.len(),
        ela.as_millis() as f64 / tests.len() as f64
    );
    println!(
        "{}\t{}\t{}",
        correct,
        test_num,
        correct as f64 / test_num as f64
    );
    println!("Naive alignments");
    println!("answer\tpredict\tNearestFrom1\tNearestfrom2");
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1: u32 = test
                .iter()
                .zip(dataset1.iter())
                .filter_map(|(test, d1)| {
                    d1.iter().map(|seq| edlib_sys::global_dist(seq, test)).min()
                })
                .sum();
            let l2: u32 = test
                .iter()
                .zip(dataset2.iter())
                .filter_map(|(test, d2)| {
                    d2.iter().map(|seq| edlib_sys::global_dist(seq, test)).min()
                })
                .sum();
            let p = if l2 < l1 { 2 } else { 1 };
            println!("{}\t{}\t{}\t{}", ans, p, l1, l2);
            *ans == p
        })
        .count();
    println!(
        "{}\t{}\t{}",
        correct,
        test_num,
        correct as f64 / test_num as f64
    );
}

fn generate_dataset<T: Rng>(
    templates: &[Vec<u8>],
    min_num: usize,
    max_num: usize,
    rng: &mut T,
    k: usize,
) -> (Vec<Vec<Vec<u8>>>, Vec<DBGHMM>) {
    let dataset: Vec<_> = templates
        .iter()
        .map(|e| {
            let num = rng.gen_range(min_num, max_num);
            (0..num)
                .map(|_| gen_sample::introduce_randomness(e, rng, &gen_sample::PROFILE))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut f = Factory::new();
    let models: Vec<_> = dataset.iter().map(|e| f.generate(e, k)).collect();
    (dataset, models)
}

fn predict(models: &[DBGHMM], test: &[Vec<u8>]) -> Vec<f64> {
    models
        .iter()
        .zip(test.iter())
        .map(|(e, f)| e.forward(f, &DEFAULT_CONFIG))
        .collect()
}

fn merge_predict(l1: &[f64], l2: &[f64]) -> (f64, f64) {
    //println!("merging prediction below:");
    let ratio: Vec<_> = l1
        .iter()
        .zip(l2.iter())
        .map(|(&x1, &x2)| as_weight(x1, x2))
        .collect();
    let weights: Vec<_> = ratio
        .iter()
        .map(|(x1, x2)| 2f64.ln() + x1 * x1.ln() + x2 * x2.ln())
        .collect();
    let tot = weights.iter().fold(0., |x, y| x + y);
    eprintln!("Dump weights");
    for (idx, ((&l1, &l2), w)) in l1.iter().zip(l2.iter()).zip(weights.iter()).enumerate() {
        let (w1, w2) = as_weight(l1, l2);
        eprintln!("{}\t{:.4}\t{:.4}\t{:.4}", idx, w1, w2, w / tot);
    }
    let (p1, p2) = ratio
        .into_iter()
        .zip(weights.into_iter())
        .map(|((f1, f2), w)| (f1 * w, f2 * w))
        .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
    assert!((p1 / tot + p2 / tot - 1.).abs() < 0.0001);
    eprintln!("Summping to:{}\t{}",p1 / tot, p2 / tot);
    eprintln!();
    (p1 / tot, p2 / tot)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}
