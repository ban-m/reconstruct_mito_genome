extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use create_simulation_data::construct_with_weights;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
#[macro_use]
extern crate log;
extern crate env_logger;
const K: usize = 6;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let chain_len = 20;
    let len = 150;
    let coverage = 20;
    let p = &gen_sample::Profile {
        sub: 0.002,
        ins: 0.002,
        del: 0.002,
    };
    let patt = 6;
    let seed = 1342374;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut templates = vec![];
    let template = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    templates.push(template.clone());
    templates.push(template.clone());
    for _ in 1..patt {
        templates.push(
            template
                .iter()
                .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                .collect::<Vec<_>>(),
        );
    }
    let dataset: Vec<_> = templates
        .iter()
        .flat_map(|seq| {
            (0..coverage)
                .map(|_| {
                    seq.iter()
                        .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let ws: Vec<_> = (0..patt)
        .map(|e| {
            (0..dataset.len())
                .map(|i| {
                    if (e * coverage <= i) && (i < (e + 1) * coverage) {
                        1.
                    } else {
                        0.
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    for i in 0..patt {
        for j in (i+1)..patt {
            let diff = templates[i]
                .iter()
                .zip(templates[j].iter())
                .map(|(x, y)| edlib_sys::global_dist(x, y))
                .sum::<u32>();
            debug!("Dist {} out of {} ({},{})", diff, chain_len, i, j);
            let (lkind, lkmerge) = calculate_lk_gain_by_merging(&dataset, i, j, &ws, chain_len);
            let lkgain = lkmerge - lkind;
            debug!("{}\t{}\t{:.4}\t{:.4}\t{:.4}", i, j, lkind, lkmerge, lkgain);
        }
    }
}

fn calculate_lk_gain_by_merging(
    data: &[Vec<Vec<u8>>],
    i: usize,
    j: usize,
    weights: &[Vec<f64>],
    len: usize,
) -> (f64, f64) {
    // Independent LK
    let tot = weights.iter().map(|e| e.iter().sum::<f64>()).sum::<f64>();
    let ws = weights
        .iter()
        .map(|e| e.iter().sum::<f64>() / tot)
        .collect::<Vec<_>>();
    let ind = calculate_likelihood(data, weights, &ws, len);
    // Merge
    let mut merged_weights = weights
        .iter()
        .enumerate()
        .filter_map(|(idx, e)| if idx != j { Some(e) } else { None })
        .map(|e| e.to_vec())
        .collect::<Vec<_>>();
    for (idx, w) in weights[j].iter().enumerate() {
        merged_weights[i][idx] += w;
    }
    let ws = merged_weights
        .iter()
        .map(|e| e.iter().sum::<f64>() / tot)
        .collect::<Vec<_>>();
    let merge = calculate_likelihood(data, &merged_weights, &ws, len);
    (ind, merge)
}

fn calculate_likelihood(
    data: &[Vec<Vec<u8>>],
    weights: &[Vec<f64>],
    ws: &[f64],
    _len: usize,
) -> f64 {
    let models: Vec<_> = weights
        .iter()
        .map(|weight| construct_with_weights(data, weight, K))
        .collect();
    data.par_iter()
        .map(|read| {
            let ms: Vec<_> = models
                .iter()
                .zip(ws.iter())
                .map(|(ms, w)| {
                    read.iter()
                        .zip(ms.iter())
                        .map(|(chunk, model)| model.forward(chunk, &DEFAULT_CONFIG))
                        .sum::<f64>()
                        + w.ln()
                })
                .collect();
            logsumexp(&ms)
        })
        .sum::<f64>()
}

// log sum_i exp xs[i]
fn logsumexp(xs: &[f64]) -> f64 {
    let max = xs
        .iter()
        .fold(std::f64::MIN, |x, &y| if x < y { y } else { x });
    max + xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln()
}
