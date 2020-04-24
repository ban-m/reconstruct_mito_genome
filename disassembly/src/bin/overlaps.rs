extern crate bio_utils;
use bio_utils::fasta;
use bio_utils::lasttab;
use std::collections::HashMap;
use std::fs;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut references: HashMap<_, _> = fasta::parse_into_vec(&args[1])?
        .into_iter()
        .map(|record| (record.id().to_string(), vec![false; record.seq().len()]))
        .collect();
    use std::io::{BufRead, BufReader};
    let alignments: Vec<_> = fs::File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| lasttab::LastTAB::from_line(&e))
        .collect();
    let threshold: f64 = (0.1f64).powi(50);
    for alignment in alignments.into_iter().filter(|l| l.e_score() < threshold) {
        //println!("{}", alignment);
        if let Some(reference) = references.get_mut(alignment.seq1_name()) {
            let start = alignment.seq1_start_from_forward();
            let end = alignment.seq1_end_from_forward();
            for i in start..end {
                reference[i] = true;
            }
        }
    }
    for (id, reference) in references.into_iter() {
        let len = reference.len() as f32;
        let cover_ratio = reference.iter().map(|&e| (e as u32) as f32).sum::<f32>() / len;
        println!("{}\t{}", id, cover_ratio);
    }
    Ok(())
}
