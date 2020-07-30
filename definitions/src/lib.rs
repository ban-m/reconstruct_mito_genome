//! Definitions -- A tiny interface for HLA-typing problem.
//! Roughly speaking, we incorporate with other programs, pass messages, or interact with other CLI via JSON object format. Specifically, the message is encoded only one, possibly large, structure named [DataSet](DataSet)

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DataSet {
    pub raw_reads: Vec<RawRead>,
    pub hic_pairs: Vec<HiCPair>,
    pub selected_chunks: Vec<Unit>,
    pub encoded_reads: Vec<EncodedRead>,
    pub hic_edges: Vec<HiCEdge>,
    pub assignments: Vec<Assignment>,
}

impl DataSet {
    pub fn with_param(
        raw_reads: Vec<RawRead>,
        hic_pairs: Vec<HiCPair>,
        selected_chunks: Vec<Unit>,
        encoded_reads: Vec<EncodedRead>,
        hic_edges: Vec<HiCEdge>,
        assignments: Vec<Assignment>,
    ) -> Self {
        Self {
            raw_reads,
            hic_pairs,
            selected_chunks,
            encoded_reads,
            hic_edges,
            assignments,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawRead {
    pub name: String,
    pub desc: String,
    pub id: u64,
    pub seq: String,
}

impl RawRead {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
}

impl std::fmt::Display for RawRead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} {} {}\n{}", self.name, self.desc, self.id, self.seq)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HiCPair {
    pub pair1: u64,
    pub pair2: u64,
    pub pair_id: u64,
    pub seq1: String,
    pub seq2: String,
}

impl HiCPair {
    pub fn seq1(&self) -> &[u8] {
        self.seq1.as_bytes()
    }
    pub fn seq2(&self) -> &[u8] {
        self.seq2.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Unit {
    pub id: u64,
    pub seq: String,
}

impl Unit {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EncodedRead {
    pub original_length: usize,
    pub leading_gap: usize,
    pub trailing_gap: usize,
    pub id: u64,
    pub edges: Vec<Edge>,
    pub nodes: Vec<Node>,
}

impl EncodedRead {
    pub fn is_gappy(&self) -> bool {
        self.nodes.is_empty()
    }
    pub fn encoded_rate(&self) -> f64 {
        let encoded_length = self.encoded_length();
        encoded_length as f64 / self.original_length as f64
    }
    pub fn encoded_length(&self) -> usize {
        let sum = self.nodes.iter().map(|n| n.query_length()).sum::<usize>();
        let offset = self
            .edges
            .iter()
            .map(|e| e.offset)
            .filter(|&e| e < 0)
            .sum::<i64>();
        let length = sum as i64 + offset;
        if length < 0 {
            0
        } else {
            length as usize
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Edge {
    pub from: u64,
    pub to: u64,
    pub offset: i64,
    pub label: String,
}

impl Edge {
    pub fn label(&self) -> &[u8] {
        self.label.as_bytes()
    }
    pub fn from_nodes(ns: &[Node], seq: &[u8]) -> Self {
        let (from, to) = match *ns {
            [ref from, ref to] => (from, to),
            _ => unreachable!(),
        };
        let end = from.position_from_start + from.query_length();
        let start = to.position_from_start;
        let label = if start < end {
            "".to_string()
        } else {
            String::from_utf8_lossy(&seq[end..start]).to_string()
        };
        Edge {
            from: from.unit,
            to: to.unit,
            offset: start as i64 - end as i64,
            label,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Node {
    /// 0-index.
    pub position_from_start: usize,
    pub unit: u64,
    pub cluster: u64,
    pub seq: String,
    pub is_forward: bool,
    pub cigar: Vec<Op>,
}

impl Node {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
    pub fn query_length(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                Op::Match(l) | Op::Ins(l) => *l,
                Op::Del(_) => 0,
            })
            .sum::<usize>()
    }
    pub fn recover(&self, unit: &Unit) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (read, unit) = (self.seq(), unit.seq());
        println!("R:{}", String::from_utf8_lossy(read));
        println!("Q:{}", String::from_utf8_lossy(unit));
        let (mut q, mut al, mut r) = (vec![], vec![], vec![]);
        let (mut q_pos, mut r_pos) = (0, 0);
        for op in self.cigar.iter() {
            match *op {
                Op::Match(l) => {
                    al.extend(
                        read[q_pos..q_pos + l]
                            .iter()
                            .zip(&unit[r_pos..r_pos + l])
                            .map(|(x, y)| if x == y { b'|' } else { b'X' }),
                    );
                    q.extend(read[q_pos..q_pos + l].iter().copied());
                    r.extend(unit[r_pos..r_pos + l].iter().copied());
                    q_pos += l;
                    r_pos += l;
                }
                Op::Del(l) => {
                    al.extend(vec![b' '; l]);
                    q.extend(vec![b' '; l]);
                    r.extend(unit[r_pos..r_pos + l].iter().copied());
                    r_pos += l;
                }
                Op::Ins(l) => {
                    al.extend(vec![b' '; l]);
                    q.extend(read[q_pos..q_pos + l].iter().copied());
                    r.extend(vec![b' '; l]);
                    q_pos += l;
                }
            }
        }
        (q, al, r)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Copy)]
pub enum Op {
    Match(usize),
    /// Deletion with respect to the reference.
    Del(usize),
    /// Insertion with respect to the reference.
    Ins(usize),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HiCEdge {
    pub pair_id: u64,
    pub pair1: u64,
    pub pair2: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assignment {
    pub id: u64,
    pub cluster: usize,
}

impl Assignment {
    pub fn new(id: u64, cluster: usize) -> Self {
        Self { id, cluster }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
