//! A crate to refresent a peak call format. It is written as TSV,
//! and each record consist of "contig name(String)", "start position(int)", "# of reads stopped", and "# of reads started".
//! Note that the end position of a unit is implicitly determined by its format.
//! For example, when we have a peak call format like below:
//! ```tsv
//! ctg0   0     1000   10   40
//! ctg0   1000  2000   20   10
//! ```
//! The unit is ctg0:[0-1000) and ctg0:[1000-2000) of the contig.
//! Note that the file should be supplied with its contigs to
//! make the encoding consistent.
use bio_utils::fasta;

use std::collections::BTreeMap;
#[derive(Deserialize, Serialize, Debug, Default)]
pub struct UnitDefinitions {
    // This is the index of the contig.
    contig_index: BTreeMap<String, u16>,
    contigs: Vec<fasta::Record>,
    units: Vec<Unit>,
}

impl UnitDefinitions {
    pub fn open_peak_with_contigs(peak_file: String, contigs: Vec<fasta::Record>) -> Self {
        let contig_index: BTreeMap<_, u16> = contigs
            .iter()
            .enumerate()
            .map(|(idx, e)| (e.id().to_string(), idx as u16))
            .collect();
        let mut units:Vec<_> = peak_file
            .lines()
            .filter_map(|e| {
                let mut e = e.split('\t');
                let idx = contig_index[e.next()?];
                let start = e.next().and_then(|e| e.parse().ok())?;
                let end = e.next().and_then(|e| e.parse().ok())?;
                Some(Unit::new(idx, start, end))
            })
            .collect();
        units.sort();
        Self {
            contig_index,
            contigs,
            units,
        }
    }
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Eq, PartialEq)]
pub struct Unit {
    contig: u16,
    start: usize,
    end: usize,
}

impl Unit {
    fn new(contig: u16, start: usize, end: usize) -> Self {
        Self { contig, start, end }
    }
}

impl std::cmp::Ord for Unit {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        match self.contig.cmp(&other.contig) {
            Equal => self.start.cmp(&other.start),
            x => x,
        }
    }
}

impl std::cmp::PartialOrd for Unit {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
