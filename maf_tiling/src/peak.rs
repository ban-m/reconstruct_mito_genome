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
/// The struct to define units.
/// Inside this struct, reference contigs and
/// definition of units are there.
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
        let mut units: Vec<_> = peak_file
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
    /// Return the reference to the fasta record with specified ID. If there's no
    /// entry for the query, return None.
    pub fn get_reference_sequence(&self, id: &str) -> Option<&fasta::Record> {
        self.contig_index
            .get(id)
            .and_then(|&i| self.contigs.get(i as usize))
    }
    /// Return the definition of most nearest unit with smaller index.
    pub fn definition(&self, id: &str, position: usize) -> Option<Unit> {
        let contig = *self.contig_index.get(id)?;
        let start = position;
        let end = start;
        let probe = Unit { contig, start, end };
        match self.units.binary_search(&probe) {
            Ok(res) => self.units.get(res).cloned(),
            Err(res) => self.units.get(res).cloned(),
        }
    }
    /// Pull the determined units. If it is not the encoded one,
    /// it returns None.
    pub fn pull_unit(&self, unit: &Unit) -> Option<&[u8]> {
        self.contigs
            .get(unit.contig() as usize)
            .map(|fasta| &fasta.seq()[unit.start()..unit.end()])
    }
}

/// The definition of the unit.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Eq, PartialEq)]
pub struct Unit {
    contig: u16,
    start: usize,
    end: usize,
}

impl Unit {
    pub fn new(contig: u16, start: usize, end: usize) -> Self {
        Self { contig, start, end }
    }
    pub fn contig(&self) -> u16 {
        self.contig
    }
    pub fn start(&self) -> usize {
        self.start
    }
    pub fn end(&self) -> usize {
        self.end
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
