extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128StarStar;
use rayon::prelude::*;
fn main() {
    let len = 50;
    let num_seq = 30;
    let test_num = 1000;
    let k = 6;
    let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(12_218_993);
    let num_templates: Vec<_> = (2..10).collect();
    println!("Entropy\tSkew\tNumTemplate");
    for num_template in num_templates {
        let template = dbg_hmm::gen_sample::generate_seq(&mut rng, len);
        let templates: Vec<_> = (0..num_template)
            .map(|i| match i % 3 {
                0 => gen_sample::introduce_errors(&template, &mut rng, 1, 0, 0),
                1 => gen_sample::introduce_errors(&template, &mut rng, 0, 1, 0),
                _ => gen_sample::introduce_errors(&template, &mut rng, 0, 0, 1),
            })
            .collect();
        let p = &gen_sample::PROFILE;
        let dataset: Vec<Vec<_>> = templates
            .iter()
            .map(|t| {
                (0..num_seq)
                    .map(|_| gen_sample::introduce_randomness(t, &mut rng, p))
                    .collect()
            })
            .collect();
        let weight = vec![1.; num_seq];
        let models: Vec<_> = dataset
            .iter()
            .map(|chunks| {
                let data: Vec<_> = chunks.iter().map(|e| e.as_slice()).collect();
                DBGHMM::new_with_weight(&data, &weight, k)
            })
            .collect();
        use rand::prelude::SliceRandom;
        let tests: Vec<_> = (0..test_num)
            .map(|_| {
                let t = templates.choose(&mut rng).unwrap();
                gen_sample::introduce_randomness(&t, &mut rng, p)
            })
            .collect();
        let entropies: Vec<_> = tests
            .par_iter()
            .map(|test| entropy(&models, test))
            .collect();
        for e in entropies {
            println!("{:.4}\tSkew\t{}", e, num_template);
        }
        let dataset: Vec<Vec<_>> = (0..num_template)
            .map(|_| {
                (0..num_seq)
                    .map(|_| gen_sample::introduce_randomness(&template, &mut rng, p))
                    .collect()
            })
            .collect();
        let models: Vec<_> = dataset
            .iter()
            .map(|chunks| {
                let data: Vec<_> = chunks.iter().map(|e| e.as_slice()).collect();
                DBGHMM::new_with_weight(&data, &weight, k)
            })
            .collect();
        let tests: Vec<_> = (0..test_num)
            .map(|_| gen_sample::introduce_randomness(&template, &mut rng, p))
            .collect();
        let entropies: Vec<_> = tests
            .par_iter()
            .map(|test| entropy(&models, test))
            .collect();
        for e in entropies {
            println!("{:.4}\tNeutral\t{}", e, num_template);
        }
    }
}

fn entropy(ms: &[DBGHMM], q: &[u8]) -> f64 {
    let lks: Vec<_> = ms.iter().map(|m| m.forward(q, &DEFAULT_CONFIG)).collect();
    let sum = last_decompose::utils::logsumexp(&lks);
    let weights: Vec<_> = lks.iter().map(|lk| (lk - sum).exp()).collect();
    last_decompose::entropy(&weights)
}
