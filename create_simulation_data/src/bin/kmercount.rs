extern crate bio_utils;
use bio_utils::fasta;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let chunk_len: usize = args[2].parse().unwrap();
    let records: Vec<_> = fasta::parse_into_vec(&args[1]).unwrap();
    let k = 6;
    for (idx, record) in records.iter().enumerate() {
        use std::collections::HashSet;
        for (jdx, chunk) in record.seq().chunks(chunk_len).enumerate() {
            let kmer: HashSet<_> = chunk.windows(k).map(|e| e.to_vec()).collect();
            println!("{}:{}:{}", idx, jdx, kmer.len());
        }
    }
}
