extern crate bio_utils;
extern crate last_tiling;
use std::collections::HashSet;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignment = last_tiling::parse_tab_file(&args[3])?;
    let reads = last_tiling::encoding(&reads, &reference, &alignment);
    let unique_units: HashSet<_> = reads
        .iter()
        .flat_map(|e| e.seq())
        .filter_map(|e| e.encode())
        .map(|e| (e.contig, e.unit))
        .collect();
    println!("{}K", unique_units.len() * last_tiling::UNIT_SIZE / 1_000);
    Ok(())
}
