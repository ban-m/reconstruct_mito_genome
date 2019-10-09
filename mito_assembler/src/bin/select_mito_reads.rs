// Filtering out the reads mapped to genomic region.
// It leave the unmapped reads as-is so you might have to filter out low quality reads afterwords.
// Current "mapped" criteria is "Mapped region is more than 80 percent of entire read."
extern crate bio_utils;
use bio_utils::fasta;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
const THR:usize = 5_000;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let stdin = std::io::stdin();
    let reads_id: HashSet<_> = BufReader::new(stdin.lock())
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(is_good_aln)
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

// Determine whether or not alignment is good(i.e. cover more than THR fraction).
fn is_good_aln(paf: String) -> Option<String> {
    let contents: Vec<&str> = paf.split('\t').collect();
    let start: usize = contents[2].parse().unwrap();
    let end: usize = contents[3].parse().unwrap();
    if contents[5] == "NC_037304.1" && (end - start) > THR{
        Some(contents[0].to_string())
    } else {
        None
    }
}
