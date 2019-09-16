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
pub const SUBUNIT_SIZE: usize = 200;
use std::collections::BTreeMap;
/// The struct to define units.
/// Inside this struct, reference contigs and
/// definition of units are there.
#[derive(Deserialize, Serialize, Debug, Default)]
pub struct UnitDefinitions {
    // This is the index of the contig.
    contig_index: BTreeMap<String, u16>,
    contigs: Vec<fasta::Record>,
    peaks: Vec<Peak>,
}

impl std::fmt::Display for UnitDefinitions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for c in &self.contigs {
            writeln!(f, "ID:{}\tLEN:{}", c.id(), c.seq().len())?;
        }
        writeln!(f, "Contig name\tStart\tEnd\tNumber")?;
        for peak in &self.peaks {
            let c = &self.contigs[peak.contig as usize];
            writeln!(f, "{}\t{}\t{}\t{}", c.id(), peak.start, peak.end, peak.num)?;
        }
        Ok(())
    }
}

impl UnitDefinitions {
    pub fn open_peak_with_contigs(peak_file: String, contigs: Vec<fasta::Record>) -> Self {
        let contig_index: BTreeMap<_, u16> = contigs
            .iter()
            .enumerate()
            .map(|(idx, e)| (e.id().to_string(), idx as u16))
            .collect();
        let mut raw_peaks: Vec<(u16, usize)> = peak_file
            .lines()
            .skip(1)
            .filter_map(|e| {
                debug!("{}", e);
                let mut e = e.split('\t');
                let idx = contig_index[e.next()?];
                let start = e.next().and_then(|e| e.parse().ok())?;
                Some((idx, start))
            })
            .collect();
        raw_peaks.sort();
        let mut contigs_peaks: BTreeMap<u16, u16> =
            contig_index.values().cloned().map(|e| (e, 0)).collect();
        debug!("{:?}", contigs_peaks);
        let mut peaks: Vec<_> = raw_peaks
            .windows(2)
            .map(|w| {
                let prev = w[0];
                let next = w[1];
                let start = prev.1;
                let end = if prev.0 == next.0 {
                    next.1
                } else {
                    contigs[prev.0 as usize].len()
                };
                let num = contigs_peaks[&prev.0];
                let p = Peak::new(prev.0, start, end, num);
                let x = contigs_peaks.get_mut(&prev.0).unwrap();
                *x += 1;
                p
            })
            .collect();
        let &(idx, start) = raw_peaks.last().unwrap();
        let end = contigs[idx as usize].len();
        peaks.push(Peak::new(idx, start, end, contigs_peaks[&idx]));
        Self {
            contig_index,
            contigs,
            peaks,
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
    pub fn definition(&self, id: &str, position: usize) -> Option<Peak> {
        let contig = *self.contig_index.get(id)?;
        let start = position;
        let end = start;
        let num = 0;
        let probe = Peak::new(contig, start, end, num);
        match self.peaks.binary_search(&probe) {
            Ok(res) => self.peaks.get(res).cloned(),
            Err(res) => self.peaks.get(res - 1).cloned(),
        }
    }
    /// Pull the determined units. If it is not the encoded one,
    /// it returns None.
    pub fn pull_unit(&self, unit: &Peak) -> Option<&[u8]> {
        self.contigs
            .get(unit.contig() as usize)
            .map(|fasta| &fasta.seq()[unit.start()..unit.end()])
    }
}

/// The definition of the unit.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Eq, PartialEq)]
pub struct Peak {
    contig: u16,
    start: usize,
    end: usize,
    num: u16,
}

impl Peak {
    pub fn new(contig: u16, start: usize, end: usize, num: u16) -> Self {
        Self {
            contig,
            num,
            end,
            start,
        }
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
    pub fn num(&self) -> u16 {
        self.num
    }
}

impl std::cmp::Ord for Peak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        match self.contig.cmp(&other.contig) {
            Equal => self.start.cmp(&other.start),
            x => x,
        }
    }
}

impl std::cmp::PartialOrd for Peak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
