//! A module to represent encoded reads.
use super::lasttab;
use std::fmt;

/// A struct to represent encoded read.
/// It should be used with the corresponding UnitDefinitions.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct EncodedRead {
    pub id: String,
    pub seq: Vec<ChunkedUnit>,
}

impl EncodedRead {
    pub fn from(id: String, seq: Vec<ChunkedUnit>) -> Self {
        Self { id, seq }
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn seq(&self) -> &[ChunkedUnit] {
        &self.seq
    }
    // pub fn to_forward(&self)->Self{
    // }
    // pub fn fill_gap(&self)->Self{
    // }
}

impl fmt::Display for EncodedRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, ">{}", self.id)?;
        for unit in &self.seq {
            write!(f, "{} ", unit)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ChunkedUnit {
    En(Encode),
    Gap(GapUnit),
}

impl ChunkedUnit {
    pub fn is_gap(&self) -> bool {
        match self {
            ChunkedUnit::En(_) => false,
            ChunkedUnit::Gap(_) => true,
        }
    }
    pub fn is_encode(&self) -> bool {
        match self {
            ChunkedUnit::En(_) => true,
            ChunkedUnit::Gap(_) => false,
        }
    }
}

impl fmt::Display for ChunkedUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::En(encode) => write!(f, "Encode({})", encode),
            Self::Gap(gap) => write!(f, "Gap({})", gap),
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct GapUnit {
    bases: String,
}

impl fmt::Display for GapUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.bases.len())
    }
}

impl fmt::Debug for GapUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", &self.bases)
    }
}

impl GapUnit {
    pub fn new(seq: &[u8]) -> Self {
        let bases = String::from_utf8(seq.to_vec()).unwrap();
        Self { bases }
    }
    pub fn len(&self) -> usize {
        self.bases.len()
    }
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
    pub fn set_bases(&mut self, seq: &[u8]) {
        self.bases.clear();
        self.bases.push_str(&String::from_utf8_lossy(seq));
    }
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Encode {
    pub contig: u16,
    pub unit: u16,
    pub bases: String,
    pub ops: Vec<lasttab::Op>,
    pub is_forward: bool,
}

impl fmt::Display for Encode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let t = if self.is_forward { 'F' } else { 'R' };
        write!(f, "{}:{}({})", self.contig, self.unit, t)
    }
}

impl Encode {
    pub fn sketch(contig: u16, unit: u16, is_forward: bool) -> Self {
        let bases = String::new();
        let ops = vec![];
        Self {
            contig,
            unit,
            bases,
            ops,
            is_forward,
        }
    }
    pub fn len(&self) -> usize {
        self.bases.len()
    }
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
    pub fn set_bases(&mut self, seq: &[u8]) {
        self.bases.clear();
        self.bases.push_str(&String::from_utf8_lossy(seq));
    }
    pub fn set_ops(&mut self, ops: &[lasttab::Op]) {
        self.ops.clear();
        self.ops.extend(ops);
    }
    pub fn is_forward(&self) -> bool {
        self.is_forward
    }
    // The reference should be consistent with the `is_forward` value.
    pub fn view(&self, refr: &[u8]) {
        let (mut r, mut q) = (0, 0);
        let bases = self.bases.as_bytes();
        let (mut rs, mut qs) = (vec![], vec![]);
        for op in &self.ops {
            match op {
                lasttab::Op::Match(l) => {
                    rs.extend(&refr[r..r + l]);
                    qs.extend(&bases[q..q + l]);
                    r += l;
                    q += l;
                }
                lasttab::Op::Seq1In(l) => {
                    rs.extend(&vec![b'-'; *l]);
                    qs.extend(&bases[q..q + l]);
                    q += l;
                }
                lasttab::Op::Seq2In(l) => {
                    rs.extend(&refr[r..r + l]);
                    qs.extend(&vec![b'-'; *l]);
                    r += l;
                }
            }
        }
        println!("{}", String::from_utf8_lossy(&rs[..100]));
        println!("{}", String::from_utf8_lossy(&qs[..100]));
        println!("{}", String::from_utf8_lossy(&rs[100..]));
        println!("{}", String::from_utf8_lossy(&qs[100..]));
    }
}
