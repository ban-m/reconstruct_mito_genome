use last_tiling::unit::Encode;
///! ERread -- a lightweight version of LastEncodedRead.
use last_tiling::EncodedRead;
use last_tiling::UNIT_SIZE;
use std::fmt;
const CLIP_THR: usize = 2000;
const MARGIN: usize = 20;
/// A simple repr for EncodedRead.
#[derive(Debug, Clone)]
pub struct ERead {
    pub id: String,
    pub seq: Vec<CUnit>,
    pub has_head_clip: bool,
    pub has_tail_clip: bool,
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

impl fmt::Display for ERead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, ">{}", self.id)?;
        for unit in &self.seq {
            write!(f, "{} ", unit)?;
        }
        Ok(())
    }
}

impl ERead {
    pub fn new_no_gapfill(er: EncodedRead) -> Self {
        use last_tiling::unit::ChunkedUnit;
        let EncodedRead { id, seq } = er;
        // Check whether it has head clip.
        let has_head_clip = match seq.first() {
            Some(ChunkedUnit::Gap(ref gap)) => gap.len() > CLIP_THR,
            _ => false,
        };
        // Check whether it has head clip.
        let has_tail_clip = match seq.last() {
            Some(ChunkedUnit::Gap(ref gap)) => gap.len() > CLIP_THR,
            _ => false,
        };
        let seq: Vec<_> = seq
            .iter()
            .filter_map(|read| read.encode())
            .map(|e| CUnit::new(e.clone()))
            .collect();
        Self {
            id,
            seq,
            has_head_clip,
            has_tail_clip,
        }
    }
    pub fn new(er: EncodedRead) -> Self {
        use last_tiling::unit::ChunkedUnit;
        let EncodedRead { id, seq } = er;
        let mut e_seq = vec![];
        // Check whether it has head clip.
        let has_head_clip = match seq.first() {
            Some(ChunkedUnit::Gap(ref gap)) if gap.len() > CLIP_THR => true,
            Some(ChunkedUnit::Gap(_)) => false,
            Some(ChunkedUnit::En(ref e)) => {
                e_seq.push(CUnit::new(e.clone()));
                false
            }
            None => false,
        };
        let e_read = seq.windows(3).filter_map(|window| {
            if window[1].is_encode() {
                let c = window[1].encode().unwrap().clone();
                return Some(CUnit::new(c));
            } else if window[0].is_encode() && window[2].is_encode() {
                // These unwraps are safe.
                let prev = window[0].encode().unwrap();
                let gap = window[1].gap().unwrap();
                let next = window[2].encode().unwrap();
                let is_same_direction = prev.is_forward() == next.is_forward();
                let is_same_contig = prev.contig == next.contig;
                let is_consective_forward = prev.is_forward() && prev.unit + 2 == next.unit;
                let is_consective_reverse = !prev.is_forward() && prev.unit == next.unit + 2;
                let is_gap_size_moderate =
                    UNIT_SIZE - MARGIN <= gap.len() && gap.len() <= UNIT_SIZE + MARGIN;
                if is_same_direction
                    && is_same_contig
                    && (is_consective_forward || is_consective_reverse)
                    && is_gap_size_moderate
                {
                    // we can fill the gap!
                    debug!("Fill gap!");
                    let contig = prev.contig;
                    let unit = if prev.is_forward() {
                        prev.unit + 1
                    } else {
                        prev.unit - 1
                    };
                    let bases = if prev.is_forward() {
                        gap.bases().to_vec()
                    } else {
                        last_tiling::revcmp(gap.bases())
                    };
                    return Some(CUnit {
                        contig,
                        unit,
                        bases,
                    });
                }
            }
            None
        });
        e_seq.extend(e_read);
        let has_tail_clip = match seq.last() {
            Some(ChunkedUnit::Gap(ref gap)) if gap.len() > CLIP_THR => true,
            Some(ChunkedUnit::Gap(_)) => false,
            Some(ChunkedUnit::En(ref e)) => {
                e_seq.push(CUnit::new(e.clone()));
                false
            }
            None => false,
        };
        Self {
            id,
            seq: e_seq,
            has_head_clip,
            has_tail_clip,
        }
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
        let (has_head_clip, has_tail_clip) = (false, false);
        Self {
            id,
            seq,
            has_head_clip,
            has_tail_clip,
        }
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
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
    pub fn has_head_clip(&self) -> bool {
        self.has_head_clip
    }
    pub fn has_tail_clip(&self) -> bool {
        self.has_tail_clip
    }
    pub fn does_touch(&self, contig: u16, start: u16, end: u16) -> bool {
        self.seq
            .iter()
            .filter(|unit| unit.contig == contig)
            .any(|unit| start < unit.unit && unit.unit < end)
    }
    pub fn seq_mut(&mut self) -> &mut Vec<CUnit> {
        self.seq.as_mut()
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
    pub fn new(e: Encode) -> Self {
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
    }
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
    pub fn len(&self) -> usize {
        self.bases.len()
    }
    pub fn unit(&self) -> usize {
        self.unit as usize
    }
    pub fn contig(&self) -> usize {
        self.contig as usize
    }
}

impl fmt::Display for CUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}", self.contig, self.unit)
    }
}
