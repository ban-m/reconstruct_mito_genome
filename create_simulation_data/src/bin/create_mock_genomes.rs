extern crate bio_utils;
extern crate poa_hmm;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::fs::File;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reference_len: usize = args[1].parse().unwrap();
    let mut_rate: f64 = args[2].parse().unwrap();
    let _number: usize = args[3].parse().unwrap();
    let seed: u64 = args.get(4).and_then(|e| e.parse().ok()).unwrap_or(213123);
    let outname: &str = &args[5];
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let p = &gen_sample::Profile {
        sub: mut_rate / 3.,
        ins: mut_rate / 3.,
        del: mut_rate / 3.,
    };
    let template: Vec<_> = (0..4)
        .map(|_| gen_sample::generate_seq(&mut rng, reference_len / 4))
        .collect();
    let record1 = {
        let desc = Some("depth=1.0 circular=true".to_string());
        let seq: Vec<_> = template.iter().flat_map(|e| e.iter()).copied().collect();
        bio_utils::fasta::Record::with_data(&format!("reference:{}", seed), &desc, &seq)
    };
    let looped = {
        let seq: Vec<_> = template[0]
            .iter()
            .chain(template[3].iter())
            .chain(template[2].iter())
            .copied()
            .collect();
        let desc = Some("depth=0.8 circular=false".to_string());
        let seq = gen_sample::introduce_randomness(&seq, &mut rng, p);
        bio_utils::fasta::Record::with_data(&format!("looped:{}", seed), &desc, &seq)
    };
    let mut wtr = fasta::Writer::new(File::create(outname.to_string() + "_reference.fa")?);
    wtr.write_record(&record1)?;
    let mut wtr = fasta::Writer::new(File::create(outname.to_string() + "_contigs.fa")?);
    wtr.write_record(&record1)?;
    wtr.write_record(&looped)
}
