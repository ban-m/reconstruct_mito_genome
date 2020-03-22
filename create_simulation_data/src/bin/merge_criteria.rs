extern crate create_simulation_data;
extern crate edlib_sys;
extern crate last_decompose;
extern crate poa_hmm;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use last_decompose::ERead;
use poa_hmm::gen_sample;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
#[macro_use]
extern crate log;
extern crate env_logger;
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let chain_len = 40;
    let len = 150;
    let (coverage, clusters, seed, ws) = if args.len() > 3 {
        let coverage: usize = args[1].parse().unwrap();
        let clusters: usize = args[2].parse().unwrap();
        let seed: u64 = args[3].parse().unwrap();
        let ws: Vec<f64> = args[4..4 + clusters]
            .iter()
            .filter_map(|e| e.parse().ok())
            .collect();
        (coverage, clusters, seed, ws)
    } else {
        (60, 6, 1342374, vec![6f64.recip(); 6])
    };
    let p = &gen_sample::Profile {
        sub: 0.00005,
        ins: 0.00005,
        del: 0.00005,
    };
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let template = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    let templates: Vec<_> = (0..clusters)
        .map(|_| {
            template
                .iter()
                .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                .collect::<Vec<_>>()
        })
        .collect();
    let (dataset, answers): (Vec<_>, Vec<_>) = templates
        .iter()
        .zip(ws.clone())
        .enumerate()
        .flat_map(|(idx, (seq, w))| {
            (0..(coverage as f64 * w).floor() as usize)
                .map(|i| {
                    let seq = seq
                        .iter()
                        .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                        .collect::<Vec<_>>();
                    let id = format!("{:03}{:03}", idx, i);
                    (ERead::new_with_lowseq(seq, &id), idx)
                })
                .collect::<Vec<_>>()
        })
        .unzip();
    // let answers: Vec<_> = (0..dataset.len()).map(|i| (i / coverage) as u8).collect();
    let weights_of_reads: Vec<_> = answers
        .iter()
        .map(|&e| {
            let mut ws = vec![0.; clusters];
            ws[e as usize] = 1.;
            ws
        })
        .collect();
    let c = &poa_hmm::DEFAULT_CONFIG;
    let wor = &weights_of_reads;
    use last_decompose::poa_clustering::DEFAULT_ALN;
    use last_decompose::poa_clustering::{likelihood_by_merging, likelihood_of_assignments};
    let before_lk = likelihood_of_assignments(&dataset, &wor, clusters, c, &DEFAULT_ALN);
    for i in 0..clusters {
        for j in (i + 1)..clusters {
            let diff = templates[i]
                .iter()
                .zip(templates[j].iter())
                .map(|(x, y)| edlib_sys::global_dist(x, y))
                .sum::<u32>();
            debug!("Dist {} out of {} ({},{})", diff, chain_len, i, j);
            let after_lk = likelihood_by_merging(&dataset, wor, i, j, clusters, c, &DEFAULT_ALN);
            let num = (ws[i] + ws[j]) * coverage as f64;
            let lkgain = (after_lk - before_lk) / num;
            debug!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                i, j, before_lk, after_lk, lkgain
            );
            println!("{}\t{}\t{}\t{}", seed, lkgain, diff, chain_len * len);
        }
    }
}
