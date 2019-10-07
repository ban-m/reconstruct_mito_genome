// Filtering out the reads mapped to genomic region.
// It leave the unmapped reads as-is so you might have to filter out low quality reads afterwords.
// Current "mapped" criteria is "Mapped region is more than 80 percent of entire read."
extern crate bio_utils;
use bio_utils::fasta;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let stdin = std::io::stdin();
    let reads_id: HashSet<_> = BufReader::new(stdin.lock())
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.split('\t').nth(0).map(|e|e.to_string()))
        .collect();
    let stdout = std::io::stdout();
    let mut stdout = fasta::Writer::new(stdout.lock());
    fasta::parse_into_vec(&args[1])?
        .into_iter()
        .filter(|e| reads_id.contains(e.id()))
        .filter_map(|read| stdout.write_record(&read).ok())
        .count();
    Ok(())
}
