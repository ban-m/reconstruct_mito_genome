//! A module to represent encoded reads.
use super::lasttab;
use std::fmt;

/// A struct to represent encoded read.
/// It should be used with the corresponding UnitDefinitions.
#[derive(Debug, Default, Clone,Serialize,Deserialize)]
pub struct EncodedRead {
    pub id: String,
    pub seq: Vec<ChunkedUnit>,
}

impl EncodedRead {
    pub fn from(id: String, seq: Vec<ChunkedUnit>) -> Self {
        Self { id, seq }
    }
    pub fn id(&self)->&str{
        &self.id
    }
    pub fn seq(&self)->&[ChunkedUnit]{
        &self.seq
    }
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

#[derive(Debug, Clone,Serialize,Deserialize)]
pub enum ChunkedUnit {
    En(Encode),
    Gap(GapUnit),
}

impl fmt::Display for ChunkedUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::En(encode) => write!(f, "Encode({})", encode),
            Self::Gap(gap) => write!(f, "Gap({})", gap),
        }
    }
}

#[derive(Clone,Serialize,Deserialize)]
pub struct GapUnit {
    bases: Vec<u8>,
}

impl fmt::Display for GapUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.bases.len())
    }
}

impl fmt::Debug for GapUnit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.bases))
    }
}

impl GapUnit {
    pub fn new(seq: &[u8]) -> Self {
        let bases = seq.to_vec();
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
        self.bases.extend(seq);
    }
}

#[derive(Debug, Default, Clone,Serialize,Deserialize)]
pub struct Encode {
    pub contig: u16,
    pub uit: u16,
    pub subunit: u16,
    pub bases: Vec<u8>,
    pub ops: Vec<lasttab::Op>,
}

impl fmt::Display for Encode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}:{}", self.contig, self.unit, self.subunit)
    }
}

impl Encode {
    pub fn sketch(contig: u16, unit: u16, subunit: u16) -> Self {
        let bases = vec![];
        let ops = vec![];
        Self {
            contig,
            unit,
            subunit,
            bases,
            ops,
        }
    }
    pub fn len(&self)->usize{
        self.bases.len()
    }
    pub fn is_empty(&self)->bool{
        self.bases.is_empty()
    }
    pub fn set_bases(&mut self, seq: &[u8]) {
        self.bases.clear();
        self.bases.extend(seq);
    }
    pub fn set_ops(&mut self, ops: &[lasttab::Op]) {
        self.ops.clear();
        self.ops.extend(ops);
    }
    pub fn view(&self, refr: &[u8]) {
        let (mut r, mut q) = (0, 0);
        let (mut rs, mut qs) = (vec![], vec![]);
        for op in &self.ops {
            match op {
                lasttab::Op::Match(l) => {
                    rs.extend(&refr[r..r + l]);
                    qs.extend(&self.bases[q..q + l]);
                    r += l;
                    q += l;
                }
                lasttab::Op::Seq1In(l) => {
                    rs.extend(&vec![b'-'; *l]);
                    qs.extend(&self.bases[q..q + l]);
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
