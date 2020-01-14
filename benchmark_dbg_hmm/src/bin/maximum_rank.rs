extern crate dbg_hmm;
extern crate rand;
extern crate rayon;
use dbg_hmm::gen_sample::*;
use rand::{rngs::StdRng, SeedableRng};
use rayon::prelude::*;
fn main() {
    let k = 6;
    // let coverage: Vec<_> = (10..200).collect();
    // let lens: Vec<_> = (50..200).collect();
    let coverage: Vec<_> = (5..100).collect();
    let lens: Vec<_> = (150..151).collect();
    let rep = 100;
    let result: Vec<_> = coverage
        .into_par_iter()
        .flat_map(|cov| {
            let mut rng: StdRng = SeedableRng::seed_from_u64(121332983 + cov as u64);
            let mut result = vec![];
            lens.iter().for_each(|&len| {
                (0..rep).for_each(|_| {
                    let template = generate_seq(&mut rng, len);
                    let mut counter: Vec<_> = vec![0; 1 << (2 * k)];
                    for seq in (0..cov).map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
                    {
                        // eprintln!("{}",String::from_utf8_lossy(&seq));
                        for kmer in seq.windows(k) {
                            counter[to_usize(kmer)] += 1;
                        }
                    }
                    let (argmin, min) = template
                        .windows(k)
                        .map(|kmer| (kmer, counter[to_usize(kmer)]))
                        .min_by_key(|e| e.1)
                        .unwrap();
                    let kmer = String::from_utf8(argmin.to_vec()).unwrap();
                    let rank = counter.iter().filter(|&&c| c > min).count();
                    result.push((len, cov, min, rank, kmer))
                })
            });
            result
        })
        .collect();
    println!("Length\tCoverage\tMin\tRank\tArgmin");
    for (len, cov, min, rank, argmin) in result {
        println!("{}\t{}\t{}\t{}\t{}", len, cov, min, rank, argmin);
    }
}

fn to_usize(kmer: &[u8]) -> usize {
    fn to_binary(&b: &u8) -> usize {
        match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        }
    }
    kmer.iter().fold(0, |x, base| (x << 2) | to_binary(base))
}
