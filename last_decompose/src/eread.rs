///! ERread -- a lightweight version of LastEncodedRead.
use last_tiling::EncodedRead;

/// A simple repr for EncodedRead.
#[derive(Debug, Clone)]
pub struct ERead {
    pub id: String,
    pub seq: Vec<CUnit>,
}

use std::hash::{Hash, Hasher};

impl Hash for ERead {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl PartialEq for ERead {
    fn eq(&self, other: &Self) -> bool {
        self.id() == other.id()
    }
}
impl Eq for ERead {}

/// Currently, gap units are ignored.
impl ERead {
    pub fn new(er: EncodedRead) -> Self {
        let EncodedRead { id, seq } = er;
        let seq = seq
            .iter()
            .filter_map(|u| u.encode())
            .map(|e| {
                let contig = e.contig;
                let unit = e.unit;
                let bases = if e.is_forward {
                    e.bases.as_bytes().to_vec()
                } else {
                    last_tiling::revcmp(e.bases.as_bytes())
                };
                CUnit {
                    contig,
                    unit,
                    bases,
                }
            })
            .collect();
        Self { id, seq }
    }
    pub fn new_with_lowseq(raw_read: Vec<Vec<u8>>, id: &str) -> Self {
        let seq = raw_read
            .into_iter()
            .enumerate()
            .map(|(unit, bases)| {
                let contig = 0;
                let unit = unit as u16;
                CUnit {
                    contig,
                    unit,
                    bases,
                }
            })
            .collect();
        let id = id.to_string();
        Self { id, seq }
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn seq(&self) -> &[CUnit] {
        &self.seq
    }
    pub fn has(&self, contig: u16) -> bool {
        self.seq.iter().any(|c| c.contig == contig)
    }
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

/// A chunk of sequence. It is "canonicalized".
/// In other words, it is reverse-complimented if needed.
#[derive(Debug, Clone)]
pub struct CUnit {
    pub contig: u16,
    pub unit: u16,
    pub bases: Vec<u8>,
}

impl CUnit {
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
    pub fn unit(&self) -> usize {
        self.unit as usize
    }
    pub fn contig(&self) -> usize {
        self.contig as usize
    }
}
