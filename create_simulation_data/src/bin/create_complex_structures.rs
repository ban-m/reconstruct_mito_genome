extern crate bio_utils;
extern crate dbg_hmm;
extern crate rand;
extern crate rand_xoshiro;
use bio_utils::fasta::Record;
use bio_utils::fasta::Writer;
use dbg_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::fs::File;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let seed = 2312789;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let p = &gen_sample::Profile {
        sub: 0.003,
        ins: 0.003,
        del: 0.003,
    };
    let reference_len: usize = args[1].parse::<usize>().unwrap() / 4;
    let outpath = &args[2];
    std::fs::create_dir_all(outpath)?;
    let seq_a = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_b = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_c = gen_sample::generate_seq(&mut rng, reference_len);
    let seq_d = gen_sample::generate_seq(&mut rng, reference_len);
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
        let mut wtr = Writer::new(File::create(&format!("{}/reference.fa",outpath))?);
        wtr.write_record(&mc)?;
    }
    let looping_structure: Vec<_> = seq_a
        .iter()
        .chain(seq_d.iter())
        .chain(seq_c.iter())
        .copied()
        .collect();
    let looping_structure = gen_sample::introduce_randomness(&looping_structure, &mut rng, p);
    let ls = Record::with_data(
        &format!("looped:{}", seed),
        &Some("depth=0.4 circular=false".to_string()),
        &looping_structure,
    );
    {
        let mut wtr = Writer::new(File::create(&format!("{}/complex1.fa",outpath))?);
        wtr.write_record(&mc)?;
        wtr.write_record(&ls)?;
    }
    let a_to_d_circle: Vec<_> = seq_a.iter().chain(seq_d.iter()).copied().collect();
    let a_to_d_circle = gen_sample::introduce_randomness(&a_to_d_circle, &mut rng, p);
    let ad = Record::with_data(
        &format!("AtoD:{}", seed),
        &Some("depth=0.3 circular=true".to_string()),
        &a_to_d_circle,
    );
    let c_to_d_circle = seq_c
        .iter()
        .chain(seq_d.iter())
        .copied()
        .collect::<Vec<_>>();
    let c_to_d_circle = gen_sample::introduce_randomness(&c_to_d_circle, &mut rng, p);
    let cd = Record::with_data(
        &format!("CtoD:{}", seed),
        &Some("depth=0.35 circular=true".to_string()),
        &c_to_d_circle,
    );
    {
        let mut wtr = Writer::new(File::create(&format!("{}/complex2.fa",outpath))?);
        wtr.write_record(&mc)?;
        wtr.write_record(&ad)?;
        wtr.write_record(&cd)?;
    }
    Ok(())
}
