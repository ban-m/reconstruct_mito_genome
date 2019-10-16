extern crate bio_utils;
extern crate last_tiling;
use bio_utils::fasta;
use last_tiling;
const THR: usize = 1_000;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let contigs: Vec<_> = fasta::parse_into_vec(&args[1])?;
    let alns: Vec<_> = last_tiling::parse_tab_file(&args[2])?
        .into_iter()
        .filter(|e| e.seq1_match_len() > THR && e.seq2_match_len() > THR)
        .collect();
    let mut result = vec![];
    while !alns.is_empty() {
        // Count repetitiveness.
        let mut count: HashMap<_, _> = contigs
            .iter()
            .map(|e| (e.id().to_string(), vec![0; e.seq().len()]))
            .collect();
        let c1 = count.get_mut(aln.seq1_name()).unwrap();
        for aln in &alns {
            for i in seq1_start_from_forward()..aln.seq1_end_from_forward() {
                c1[i] += 1;
            }
        }
        let c2 = count.get_mut(aln.seq2_name()).unwrap();
        for aln in &alns {
            for i in seq1_start_from_forward()..aln.seq1_end_from_forward() {
                c1[i] += 1;
            }
        }
    }
}
