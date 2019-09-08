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
use std::path::Path;

pub fn parse_tab_file<P: AsRef<Path>>(tab_file: P) -> std::io::Result<Vec<LastTAB>> {
    let lines = std::fs::read_to_string(tab_file)?;
    Ok(lines.lines().filter_map(LastTAB::from_line).collect())
}

pub fn parse_peak_file<P: AsRef<Path>>(peak_file: P, ctg_file: P) -> std::io::Result<UnitDefinitions> {
    let ctgs:Vec<_> = bio_utils::fasta::parse_into_vec(ctg_file)?;
    let peaks = std::fs::read_to_string(peak_file)?;
    Ok(peak::UnitDefinitions::open_peak_with_contigs(peaks,ctgs))
}

pub fn encoding(
    fasta: &[fasta::Record],
    defs: &UnitDefinitions,
    alns: &[LastTAB],
) -> Vec<EncodedRead> {
    vec![]
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
