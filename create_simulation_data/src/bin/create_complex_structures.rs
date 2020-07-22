extern crate bio_utils;
extern crate poa_hmm;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta::Record;
use bio_utils::fasta::Writer;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::fs::File;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let seed = 2312789;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut_rate: f64 = args[3].parse().unwrap();
    let p = &gen_sample::Profile {
        sub: mut_rate / 3.,
        ins: mut_rate / 3.,
        del: mut_rate / 3.,
    };
    let reference_len: usize = args[1].parse::<usize>().unwrap() / 5;
    let outpath = &args[2];
    std::fs::create_dir_all(outpath)?;
    let seq_a = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_b = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_c = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_d = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_e = gen_sample::generate_seq(&mut rng, reference_len);
    let master_circle: Vec<_> = seq_a
        .iter()
        .chain(seq_b.iter())
        .chain(seq_c.iter())
        .chain(seq_d.iter())
        .copied()
        .collect();
    let mc = Record::with_data(
        &format!("master_circle:{}", seed),
        &Some("depth=1.0 circular=true".to_string()),
        &master_circle,
    );

    {
        let mut wtr = Writer::new(File::create(&format!("{}/reference.fa", outpath))?);
        wtr.write_record(&mc)?;
    }
    let sv_1: Vec<_> = seq_a
        .iter()
        .chain(seq_b.iter())
        .chain(seq_e.iter())
        .chain(seq_e.iter())
        .copied()
        .collect();
    let sv_1 = gen_sample::introduce_randomness(&sv_1, &mut rng, p);
    let sv_1 = Record::with_data(
        &format!("ABED:{}", seed),
        &Some("depth=0.6 circular=false".to_string()),
        &sv_1,
    );
    fn revcmp(xs: &[u8]) -> Vec<u8> {
        xs.iter()
            .rev()
            .map(|e| match e {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'T',
                _ => panic!(),
            })
            .collect::<Vec<_>>()
    };
    let sv_2 = seq_b.clone();
    let sv_2 = gen_sample::introduce_randomness(&sv_2, &mut rng, p);
    let sv_2 = Record::with_data(
        &format!("B"),
        &Some("depth=0.5 circulra=false".to_string()),
        &sv_2,
    );
    let sv_3: Vec<_> = seq_a.iter().chain(revcmp(&seq_c).iter()).copied().collect();
    let sv_3 = gen_sample::introduce_randomness(&sv_3, &mut rng, p);
    let sv_3 = Record::with_data(
        &format!("ArevC"),
        &Some("depth=0.8 cirular=false".to_string()),
        &sv_3,
    );
    {
        let mut wtr = Writer::new(File::create(&format!("{}/complex.fa", outpath))?);
        wtr.write_record(&mc)?;
        wtr.write_record(&sv_1)?;
        wtr.write_record(&sv_2)?;
        wtr.write_record(&sv_3)?;
    }
    Ok(())
}
