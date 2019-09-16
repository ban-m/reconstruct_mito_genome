//! This is a libray to tiling a fasta file into
//! chunks of contigs(units).
extern crate rmp_serde;
#[macro_use]
extern crate serde;
extern crate bio_utils;
use bio_utils::fasta;
mod peak;
pub use peak::UnitDefinitions;
mod unit;
pub use unit::EncodedRead;
mod lasttab;
use lasttab::LastTAB;
use std::collections::HashMap;
use std::path::Path;
pub fn parse_tab_file<P: AsRef<Path>>(tab_file: P) -> std::io::Result<Vec<LastTAB>> {
    let lines = std::fs::read_to_string(tab_file)?;
    Ok(lines.lines().filter_map(LastTAB::from_line).collect())
}

pub fn parse_peak_file<P: AsRef<Path>>(
    peak_file: P,
    ctg_file: P,
) -> std::io::Result<UnitDefinitions> {
    let ctgs: Vec<_> = bio_utils::fasta::parse_into_vec(ctg_file)?;
    let peaks = std::fs::read_to_string(peak_file)?;
    Ok(peak::UnitDefinitions::open_peak_with_contigs(peaks, ctgs))
}

pub fn encoding(
    fasta: &[fasta::Record],
    defs: &UnitDefinitions,
    alns: &[LastTAB],
) -> Vec<EncodedRead> {
    // Distribute alignments to each reads.
    // bucket[i] is the alignment for fasta[i].
    let buckets = distribute(fasta, alns);
    buckets
        .into_iter()
        .zip(fasta.iter())
        .map(|(bucket, seq)| into_encoding(bucket, seq, defs))
        .collect()
}

fn into_encoding(
    mut bucket:Vec<&LastTAB>,
    seq: &fasta::Record,
    defs: &UnitDefinitions,
) -> EncodedRead {
    bucket.sort_by_key(|aln|aln.seq2_start_from_forward());
    EncodedRead::default()
}

fn distribute<'a>(fasta: &[fasta::Record], alns: &'a [LastTAB]) -> Vec<Vec<&'a LastTAB>> {
    let mut alignments_bucket: Vec<Vec<&LastTAB>> = Vec::with_capacity(fasta.len());
    let id_to_idx: HashMap<_, _> = fasta
        .iter()
        .map(|e| e.id())
        .enumerate()
        .map(|(idx, id)| (id, idx))
        .collect();
    for aln in alns {
        alignments_bucket[id_to_idx[aln.seq2_name()]].push(aln);
    }
    alignments_bucket
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
