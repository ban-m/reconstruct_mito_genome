extern crate bio_utils;
extern crate poa_hmm;
extern crate rand;
extern crate rand_xoshiro;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reference_len: usize = args[1].parse().unwrap();
    let mut_rate: f64 = args[2].parse().unwrap();
    let number: usize = args[3].parse().unwrap();
    let seed: u64 = args.get(4).and_then(|e| e.parse().ok()).unwrap_or(213123);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let p = &gen_sample::Profile {
        sub: mut_rate / 3.,
        ins: mut_rate / 3.,
        del: mut_rate / 3.,
    };
    let desc = Some("depth=1.0 circular=true".to_string());
    let template1 = gen_sample::generate_seq(&mut rng, reference_len);
    let record1 =
        bio_utils::fasta::Record::with_data(&format!("sample0:{}", seed), &desc, &template1);
    let records: Vec<_> = (1..number)
        .map(|i| {
            let seq = gen_sample::introduce_randomness(&template1, &mut rng, p);
            bio_utils::fasta::Record::with_data(&format!("sample{}:{}", i, seed), &desc, &seq)
        })
        .collect();
    let stdout = std::io::stdout();
    let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
    wtr.write_record(&record1)?;
    for record in records {
        wtr.write_record(&record)?;
    }
    Ok(())
}
