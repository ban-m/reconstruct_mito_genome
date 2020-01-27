extern crate bio_utils;
extern crate dbg_hmm;
extern crate rand;
extern crate rand_xoshiro;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let seed = 2312789;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let p = &gen_sample::Profile {
        sub: 0.003,
        ins: 0.003,
        del: 0.003,
    };
    let reference_len: usize = args[1].parse().unwrap();
    let template1 = gen_sample::generate_seq(&mut rng, reference_len);
    let template2 = gen_sample::introduce_randomness(&template1, &mut rng, p);
    let stdout = std::io::stdout();
    let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
    let desc1 = Some("depth=0.9 circular=true".to_string());
    let desc2 = Some("depth=0.1 circular=true".to_string());
    let record1 =
        bio_utils::fasta::Record::with_data(&format!("sample1:{}", seed), &desc1, &template1);
    let record2 =
        bio_utils::fasta::Record::with_data(&format!("sample2:{}", seed), &desc2, &template2);
    wtr.write_record(&record1)?;
    wtr.write_record(&record2)
}
