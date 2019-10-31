extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
use dbg_hmm::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
fn main() {
    let seed = 10034234;
    let chain_len = 20;
    let k = 5;
    let len = 150;
    let max_num = 30;
    let min_num = 15;
    let test_num = 100;
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
    let p = gen_sample::Profile {
        sub: 0.005,
        ins: 0.005,
        del: 0.005,
    };
    let templates1: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates2: Vec<_> = templates1
        .iter()
        .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        .collect();
    for (idx, (t1, t2)) in templates1.iter().zip(templates2.iter()).enumerate() {
        println!("{}\t{}", idx, edlib_sys::global_dist(t1, t2));
    }
    let _template1: Vec<_> = templates1.iter().flat_map(|e| e).copied().collect();
    let _template2: Vec<_> = templates1.iter().flat_map(|e| e).copied().collect();
    let (dataset1, model1) = generate_dataset(&templates1, min_num, max_num, &mut rng, k);
    let (dataset2, model2) = generate_dataset(&templates2, min_num, max_num, &mut rng, k);
    println!("Generated Models");
    let tests: Vec<(_, Vec<_>)> = (0..test_num)
        .map(|_| {
            if rng.gen_bool(0.5) {
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
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let l1 = predict(&model1, test);
            let l2 = predict(&model2, test);
            let (l1, l2) = merge_predict(&l1, &l2);
            println!("{}\t{:.4}\t{:.4}", ans, l1, l2);
            ((l2 < l1 && *ans == 1) || (l1 < l2 && *ans == 2))
        })
        .count();
    println!(
        "{}\t{}\t{}",
        correct,
        test_num,
        correct as f64 / test_num as f64
    );
    println!("Naive alignments");
    let correct = tests
        .iter()
        .filter(|&(ans, ref test)| {
            let test: Vec<_> = test.iter().flat_map(|e| e).copied().collect();
            let l1 = dataset1
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, &test))
                .min()
                .unwrap();
            let l2 = dataset2
                .iter()
                .map(|seq| edlib_sys::global_dist(seq, &test))
                .min()
                .unwrap();
            println!("{}\t{}\t{}", ans, l1, l2);
            ((l1 < l2 && *ans == 1) || (l2 < l2 && *ans == 2))
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
) -> (Vec<Vec<u8>>, Vec<DBGHMM>) {
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
    let dataset: Vec<_> = dataset
        .into_iter()
        .map(|e| e.into_iter().flat_map(|e| e).collect())
        .collect();
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
    // for (idx, ((&l1, &l2), w)) in l1.iter().zip(l2.iter()).zip(weights.iter()).enumerate() {
    //     let (w1, w2) = as_weight(l1, l2);
    //     println!("{}\t{:.4}\t{:.4}\t{:.4}", idx, w1, w2, w / tot);
    // }
    let (p1, p2) = ratio
        .into_iter()
        .zip(weights.into_iter())
        .map(|((f1, f2), w)| (f1 * w, f2 * w))
        .fold((0., 0.), |(x, y), (a, b)| (x + a, y + b));
    assert!((p1 / tot + p2 / tot - 1.).abs() < 0.0001);
    (p1 / tot, p2 / tot)
}

fn as_weight(x1: f64, x2: f64) -> (f64, f64) {
    let max = x1.max(x2);
    let log_denominator = max + ((x1 - max).exp() + (x2 - max).exp()).ln();
    ((x1 - log_denominator).exp(), (x2 - log_denominator).exp())
}
