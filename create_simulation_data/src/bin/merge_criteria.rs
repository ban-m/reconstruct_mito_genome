extern crate create_simulation_data;
extern crate dbg_hmm;
extern crate edlib_sys;
extern crate last_decompose;
extern crate rand;
extern crate rand_xoshiro;
extern crate rayon;
use dbg_hmm::gen_sample;
use last_decompose::ERead;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
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
    let args: Vec<_> = std::env::args().collect();
    let chain_len = 20;
    let len = 150;
    let (coverage, clusters, seed) = if args.len() > 3 {
        let coverage: usize = args[1].parse().unwrap();
        let clusters: usize = args[2].parse().unwrap();
        let seed: u64 = args[3].parse().unwrap();
        (coverage, clusters, seed)
    } else {
        (100, 6, 1342374)
    };
    let p = &gen_sample::Profile {
        sub: 0.001,
        ins: 0.001,
        del: 0.001,
    };
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut templates = vec![];
    let template = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    templates.push(template.clone());
    templates.push(template.clone());
    for _ in 2..clusters {
        templates.push(
            template
                .iter()
                .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                .collect::<Vec<_>>(),
        );
    }
    let dataset: Vec<_> = templates
        .iter()
        .enumerate()
        .flat_map(|(idx, seq)| {
            (0..coverage)
                .map(|i| {
                    let seq = seq
                        .iter()
                        .map(|e| gen_sample::introduce_randomness(e, &mut rng, p))
                        .collect::<Vec<_>>();
                    let id = format!("{:03}{:03}", idx, i);
                    ERead::new_with_lowseq(seq, &id)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let answers: Vec<_> = (0..dataset.len()).map(|i| (i / coverage) as u8).collect();
    let contigs = vec![chain_len];
    let before_lk =
        last_decompose::likelihood_of_assignments(&dataset, &answers, K, clusters, &contigs);
    let gammas: Vec<_> = (0..dataset.len())
        .map(|idx| {
            let mut gamma = vec![0.; clusters];
            let cl = idx / coverage;
            gamma[cl] = 1.;
            gamma
        })
        .collect();
    for i in 0..clusters {
        for j in (i + 1)..clusters {
            let diff = templates[i]
                .iter()
                .zip(templates[j].iter())
                .map(|(x, y)| edlib_sys::global_dist(x, y))
                .sum::<u32>();
            debug!("Dist {} out of {} ({},{})", diff, chain_len, i, j);
            let after_lk = last_decompose::likelihood_by_merging(
                &dataset, &gammas, i, j, clusters, K, &contigs,
            );
            let lkgain = after_lk - before_lk;
            debug!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                i, j, before_lk, after_lk, lkgain
            );
            println!("{}\t{}\t{}", seed, lkgain, diff);
        }
    }
}
