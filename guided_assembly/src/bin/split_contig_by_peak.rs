extern crate bio;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let fasta: HashMap<_, _> = fasta::Reader::from_file(&Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| (e.id().to_string(), e))
        .collect();
    let mut splits = HashMap::new();
    let mut wtr = bio::io::fasta::Writer::new(std::io::stdout());
    for (contig, position) in BufReader::new(File::open(&Path::new(&args[2]))?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let con: Vec<_> = e.split('\t').collect();
            let contig = con[0].to_string();
            let position = con[1].parse().ok()?;
            Some((contig, position))
        })
    {
        let length = fasta[&contig].seq().len();
        let bucket = splits.entry(contig).or_insert(vec![0, length]);
        bucket.push(position);
    }
    let mut id = 0;
    for (contig, positions) in splits.iter_mut() {
        positions.sort();
        for window in positions.windows(2) {
            if window[1] == window[0] {
                continue;
            }
            let record = format!("{}", id);
            eprintln!("{}_{}_{},{}", fasta[contig].id(), window[0], window[1], id);
            wtr.write(&record, None, &fasta[contig].seq()[window[0]..window[1]])?;
            id += 1;
        }
    }
    Ok(())
}
