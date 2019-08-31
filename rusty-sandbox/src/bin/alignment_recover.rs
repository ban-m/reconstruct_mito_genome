extern crate bio;
extern crate bio_utils;
extern crate rand;
extern crate rust_htslib;
extern crate rusty_sandbox;
// use rand::{seq::SliceRandom,thread_rng};
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let sam_records: HashMap<_, _> = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    let fastq_records: HashMap<_, _> = rusty_sandbox::open_fastq_into_hashmap(&args[2])?;
    // let mut rng = thread_rng();
    let sam_record = &sam_records
        .values()
        .filter(|e| !e.is_empty())
        .nth(0)
        .unwrap()[0];
    let seq1 = &fastq_records[sam_record.q_name()];
    let seq1 = if sam_record.is_template() {
        seq1.seq().to_vec()
    } else {
        bio::alphabets::dna::revcomp(seq1.seq())
    };
    let seq2 = fastq_records[sam_record.r_name()].seq();
    let (seq1_p, ops, seq2_p) =
        bio_utils::sam::recover_alignment(&sam_record.cigar(), &seq1, seq2, sam_record.pos());
    let window = 200;
    for ((line1, line2), line3) in seq1_p
        .chunks(window)
        .zip(ops.chunks(window))
        .zip(seq2_p.chunks(window))
    {
        println!(
            "{}\n{}\n{}\n",
            String::from_utf8_lossy(line1),
            String::from_utf8_lossy(line2),
            String::from_utf8_lossy(line3)
        );
    }
    Ok(())
}
