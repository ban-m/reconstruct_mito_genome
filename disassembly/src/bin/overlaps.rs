extern crate bio_utils;
use bio_utils::fasta;
use bio_utils::lasttab;
use std::collections::HashMap;
use std::fs;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::io::{BufRead, BufReader};
    let alignments: Vec<_> = fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| lasttab::LastTAB::from_line(&e))
        .collect();
    let threshold: f64 = (0.1f64).powi(50);
    let mut references: HashMap<_, _> = fasta::parse_into_vec(&args[2])?
        .into_iter()
        .map(|record| (record.id().to_string(), vec![false; record.seq().len()]))
        .collect();
    for alignment in alignments.into_iter().filter(|l| l.e_score() < threshold) {
        if let Some(reference) = references.get_mut(alignment.seq1_name()) {
            let start = alignment.seq1_start_from_forward();
            let end = alignment.seq1_end_from_forward();
            for i in start..end {
                reference[i] = true;
            }
        }
    }
    let (tot, covered) = references
        .values()
        .map(|seq| {
            let tot = seq.len();
            let covered = seq.iter().filter(|&&b| b).count();
            (tot, covered)
        })
        .fold((0, 0), |(x, y), (a, b)| (x + a, y + b));
    print!("{}", covered as f64 / tot as f64);
    Ok(())
}
