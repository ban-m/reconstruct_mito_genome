#[macro_use]
extern crate log;
extern crate edlib_sys;
extern crate env_logger;
extern crate last_decompose;
extern crate poa_hmm;
extern crate rand;
extern crate rayon;
use poa_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let seed: u64 = std::env::args()
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap_or(10);
    //let mut rng: StdRng = SeedableRng::seed_from_u64(1219900);
    let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
    let chain_len = 40;
    let template: Vec<Vec<u8>> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let num_cluster = 2;
    let p = gen_sample::Profile {
        sub: 0.001 / 6.,
        ins: 0.001 / 6.,
        del: 0.001 / 6.,
    };
    let templates: Vec<Vec<_>> = (0..num_cluster)
        .map(|_| {
            template
                .iter()
                .map(|t| gen_sample::introduce_randomness(t, &mut rng, &p))
                .collect()
        })
        .collect();
    let dists: Vec<_> = (0..chain_len)
        .map(|pos| {
            let mut res = 0;
            for i in 0..num_cluster {
                for j in (i + 1)..num_cluster {
                    let t1 = &templates[i][pos];
                    let t2 = &templates[j][pos];
                    res += edlib_sys::global_dist(t1, t2);
                }
            }
            res
        })
        .collect();
    let coverages = vec![25; num_cluster];
    let mut gen = |ts: &Vec<Vec<u8>>| {
        ts.iter()
            .map(|t| gen_sample::introduce_randomness(&t, &mut rng, &gen_sample::PROFILE))
            .collect::<Vec<_>>()
    };
    let dataset: Vec<_> = coverages
        .iter()
        .enumerate()
        .flat_map(|(idx, &cov)| (0..cov).map(|_| gen(&templates[idx])).collect::<Vec<_>>())
        .collect();
    let score = |x, y| if x == y { 3 } else { -4 };
    let models: Vec<Vec<_>> = {
        let mut chunks = vec![vec![]; chain_len];
        for seq in dataset.iter() {
            for (idx, chunk) in seq.iter().enumerate() {
                chunks[idx].push(chunk.as_slice());
            }
        }
        let total_coverage = coverages.iter().sum::<usize>();
        let weights_of_reads = last_decompose::construct_initial_weights(
            &vec![],
            &vec![vec![]; total_coverage],
            num_cluster,
            total_coverage,
            10392032,
        );
        // let weights_of_reads: Vec<_> = coverages
        //     .iter()
        //     .enumerate()
        //     .flat_map(|(idx, &cov)| {
        //         let mut weight = vec![0.; num_cluster];
        //         weight[idx] = 1.;
        //         vec![weight.clone(); cov]
        //     })
        //     .collect();
        let mut weights = vec![vec![]; num_cluster];
        for weights_each_cluster in weights_of_reads {
            for (cluster, &w) in weights_each_cluster.iter().enumerate() {
                weights[cluster].push(w);
            }
        }
        let total_cov = coverages.iter().sum::<usize>();
        weights
            .iter()
            .map(|ws| {
                chunks
                    .iter()
                    .map(|cs| POA::generate_w_param(cs, ws, -6, -6, &score))
                    .collect()
            })
            .collect()
    };
    use last_decompose::variant_calling;
    let dataset: Vec<Vec<_>> = dataset
        .into_iter()
        .map(|read| read.into_iter().enumerate().collect())
        .collect();
    let variants = variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, true);
    for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
        println!("{}\t{}\t{:.4}\tCentrize", idx, d, w.abs());
    }
    let variants = variant_calling::variant_call_poa(&models, &dataset, &DEFAULT_CONFIG, false);
    for (idx, (d, w)) in dists.iter().zip(variants.iter()).enumerate() {
        println!("{}\t{}\t{:.4}\tNormal", idx, d, w.abs());
    }
    // for pos in 0..chain_len {
    //     for (idx, read) in test.iter().enumerate() {
    //         let &(_, ref unit) = read.iter().filter(|&&(p, _)| p == pos).nth(0).unwrap();
    //         let lks: Vec<_> = models
    //             .iter()
    //             .map(|ms| ms[pos].forward(unit, &DEFAULT_CONFIG))
    //             .map(|lk| format!("{}", lk))
    //             .collect();
    //         debug!("DUMP\t{}\t{}\t{}", idx, pos, lks.join("\t"));
    //     }
    // }
}
