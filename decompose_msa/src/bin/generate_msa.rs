extern crate rand;
use rand::prelude::*;
const BASE: &[u8; 4] = b"ACGT";
use std::io::{BufWriter, Write};
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let column_size: usize = args[1].parse().unwrap();
    let row_size: usize = args[2].parse().unwrap();
    let variant_rate: f64 = args[3].parse().unwrap();
    let mix_rate: f64 = args[4].parse().unwrap();
    let subst_rate: f64 = args[5].parse().unwrap();
    let (result, mut variant_site) =
        generate(column_size, row_size, variant_rate, mix_rate, subst_rate);
    for line in result {
        println!("{}", String::from_utf8_lossy(&line));
    }
    let mut wtr = BufWriter::new(std::fs::File::create(&args[6]).unwrap());
    variant_site.sort();
    for x in &variant_site {
        writeln!(&mut wtr, "{}", x).unwrap();
    }
}
fn generate(
    length_of_read: usize,
    number_of_read: usize,
    variation_rate: f64,
    mixing_rate: f64,
    substitution_rate: f64,
) -> (Vec<Vec<u8>>, Vec<usize>) {
    // Generate template.
    let mut rng: StdRng = SeedableRng::seed_from_u64(2328904);
    let template: Vec<u8> = (0..length_of_read)
        .filter_map(|_| BASE.choose(&mut rng))
        .copied()
        .collect();
    assert_eq!(template.len(), length_of_read);
    // Make variant.
    let variant_num = (length_of_read as f64 * variation_rate).floor() as usize;
    let variant_position: Vec<usize> = {
        let seq: Vec<_> = (0..length_of_read).collect();
        seq.choose_multiple(&mut rng, variant_num)
            .copied()
            .collect()
    };
    let mut variant = template.clone();
    for &p in &variant_position {
        variant[p] = mutate(variant[p], &mut rng);
    }
    // Construct MSA
    let mut msa = (0..number_of_read).fold(vec![], |mut res, _| {
        let is_template = rng.gen_bool(mixing_rate);
        if is_template {
            res.push(template.clone());
        } else {
            res.push(variant.clone());
        }
        res
    });
    // Introduce errors.
    for read in msa.iter_mut() {
        for base in read.iter_mut() {
            if rng.gen_bool(substitution_rate) {
                *base = mutate(*base, &mut rng);
            }
        }
    }
    (msa, variant_position)
}

fn mutate<T: RngCore>(b: u8, rng: &mut T) -> u8 {
    let b = match b {
        b'A' => b"CGT".choose(rng),
        b'C' => b"AGT".choose(rng),
        b'G' => b"ACT".choose(rng),
        b'T' => b"ACG".choose(rng),
        _ => unreachable!(),
    }
    .unwrap();
    *b
}
