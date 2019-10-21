// Filtering out the reads mapped to genomic region.
// It leave the unmapped reads as-is so you might have to filter out low quality reads afterwords.
// Current "mapped" criteria is "Mapped region is more than 80 percent of entire read."
extern crate bio_utils;
use bio_utils::fasta;
use std::collections::{HashSet,HashMap};
use std::io::{BufRead, BufReader};
const THR: f64 = 0.8;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = fasta::parse_into_vec(&args[1])?;
    let mut mapped_region: HashMap<_, _> = reads
        .iter()
        .map(|e| (e.id().to_string(), vec![false; e.seq().len()]))
        .collect();
    let stdin = std::io::stdin();
    for paf in BufReader::new(stdin.lock()).lines().filter_map(|e| e.ok()) {
        let (name, start, end) = region(&paf);
        let rs = mapped_region.get_mut(name).unwrap();
        for i in start..end {
            rs[i] = true;
        }
    }
    let contained_id: HashSet<_> = mapped_region
        .into_iter()
        .filter(|(_id, rs)| {
            let len = rs.len();
            let covered = rs.iter().filter(|&&e| e).count();
            covered as f64 / len as f64 > THR
        })
        .map(|(id, _)| id)
        .collect();
    let stdout = std::io::stdout();
    let mut stdout = fasta::Writer::new(stdout.lock());
    reads
        .into_iter()
        .filter(|e| contained_id.contains(e.id()))
        .filter_map(|read| stdout.write_record(&read).ok())
        .count();
    Ok(())
}

fn region(paf: &str) -> (&str, usize, usize) {
    let contents: Vec<&str> = paf.split('\t').collect();
    let name: &str = &contents[0];
    let start: usize = contents[2].parse().unwrap();
    let end: usize = contents[3].parse().unwrap();
    (name, start, end)
}
