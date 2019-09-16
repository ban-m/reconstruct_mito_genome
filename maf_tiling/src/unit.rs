//! A module to represent encoded reads.

use super::lasttab;

/// A struct to represent encoded read.
/// It should be used with the corresponding UnitDefinitions.
#[derive(Debug, Default, Clone)]
pub struct EncodedRead {
    id: String,
    seq: Vec<Encode>,
}

#[derive(Debug, Default, Clone)]
pub struct Encode {
    contig: u8,
    unit: u8,
    subunit: u16,
    bases: Vec<u8>,
    cigar: Vec<lasttab::Op>,
}
