#[macro_use]
extern crate serde;
extern crate env_logger;
extern crate last_decompose;
extern crate log;
extern crate serde_json;
use last_decompose::{find_breakpoint::ReadClassify, CriticalRegion};
use last_tiling::unit::ChunkedUnit;
use last_tiling::{Contigs, EncodedRead};
use std::collections::HashMap;
pub fn dump_viewer(
    results: &HashMap<String, u8>,
    reads: &[EncodedRead],
    crs: &[CriticalRegion],
    contigs: &Contigs,
) -> std::io::Result<String> {
    let clusters = summarize_clusters(contigs, reads, crs, results);
    let contigs = summarize_contig(contigs, reads);
    let reads = summarize_reads(reads, results);
    let summary = Summary {
        contigs,
        reads,
        critical_regions: crs.to_vec(),
        clusters,
    };
    Ok(serde_json::ser::to_string(&summary).unwrap())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Summary {
    contigs: Vec<Contig>,
    reads: Vec<Read>,
    critical_regions: Vec<CriticalRegion>,
    clusters: Vec<Cluster>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Contig {
    name: String,
    id: u16,
    length: usize,
    coverages: Vec<u32>,
    start_stop: Vec<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Cluster {
    contigs: Vec<Contig>,
    critical_regions: Vec<usize>,
}

fn summarize_clusters(
    contigs: &Contigs,
    reads: &[EncodedRead],
    crs: &[CriticalRegion],
    results: &HashMap<String, u8>,
) -> Vec<Cluster> {
    let mut clusters: Vec<_> = results.values().copied().collect();
    clusters.sort();
    clusters.dedup();
    clusters
        .iter()
        .map(|&cl| summarize_one_cluster(contigs, reads, crs, results, cl))
        .collect()
}

fn summarize_one_cluster(
    contigs: &Contigs,
    reads: &[EncodedRead],
    crs: &[CriticalRegion],
    results: &HashMap<String, u8>,
    cl: u8,
) -> Cluster {
    let tot = reads.len();
    let tot_len = contigs
        .get_last_units()
        .iter()
        .map(|&e| e as usize)
        .sum::<usize>();
    let thr = tot / tot_len * 2;
    let reads: Vec<_> = reads
        .iter()
        .filter(|r| results.get(r.id()).map(|&asn| asn == cl).unwrap_or(false))
        .cloned()
        .collect();
    let contigs = summarize_contig(contigs, &reads);
    let crs = crs
        .iter()
        .enumerate()
        .filter_map(|(idx, cr)| {
            let counts = reads
                .iter()
                .map(|r| last_decompose::ERead::new(r.clone()))
                .filter(|r| cr.along_with(r))
                .count();
            if counts > thr {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    Cluster {
        contigs,
        critical_regions: crs,
    }
}

fn summarize_contig(contigs: &Contigs, reads: &[EncodedRead]) -> Vec<Contig> {
    let mut cs: Vec<_> = contigs
        .names()
        .iter()
        .enumerate()
        .map(|(id, name)| {
            let id = id as u16;
            let length = contigs.get(name).unwrap().len();
            let coverages = vec![0; length / last_tiling::UNIT_SIZE + 1];
            let start_stop = vec![0; length / last_tiling::UNIT_SIZE + 1];
            let name = name.to_string();
            Contig {
                id,
                length,
                coverages,
                name,
                start_stop,
            }
        })
        .collect();
    for read in reads {
        let mut first = true;
        for unit in &read.seq {
            match unit {
                ChunkedUnit::En(encode) => {
                    if first {
                        cs[encode.contig as usize].start_stop[encode.unit as usize] += 1;
                        first = false;
                    }
                    cs[encode.contig as usize].coverages[encode.unit as usize] += 1
                }
                _ => {}
            }
        }
        if let Some(ChunkedUnit::En(encode)) =
            &read.seq.iter().rev().filter(|e| e.is_encode()).nth(0)
        {
            cs[encode.contig as usize].start_stop[encode.unit as usize] += 1;
        }
    }
    cs
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Read {
    name: String,
    units: Vec<Unit>,
    cluster: i8,
}

impl Read {
    fn from(read: &EncodedRead, assign: &HashMap<String, u8>) -> Self {
        let name = String::from(read.id());
        let cluster = assign.get(&name).map(|&e| e as i8).unwrap_or(-1);
        let mut units = vec![];
        let mut start = std::u16::MAX;
        let (mut prev_c, mut prev_u) = (std::u16::MAX, std::u16::MAX);
        for unit in read.seq().into_iter() {
            if unit.is_encode() {
                let u = unit.encode().unwrap();
                let diff = if u.unit < prev_u {
                    prev_u - u.unit
                } else {
                    u.unit - prev_u
                };
                if u.contig == prev_c && diff < 2 {
                    prev_u = u.unit;
                } else {
                    //If not right after a gap.
                    if prev_c != std::u16::MAX {
                        let e = Unit::E(prev_c, start, prev_c, prev_u);
                        units.push(e);
                    }
                    start = u.unit;
                    prev_u = u.unit;
                    prev_c = u.contig;
                }
            } else {
                // If not the head.
                if prev_c != std::u16::MAX {
                    let e = Unit::E(prev_c, start, prev_c, prev_u);
                    units.push(e);
                }
                prev_c = std::u16::MAX;
                prev_u = std::u16::MAX;
                units.push(Unit::G(unit.gap().unwrap().len()));
            }
        }
        if prev_c != std::u16::MAX && prev_u != std::u16::MAX {
            let e = Unit::E(prev_c, start, prev_c, prev_u);
            units.push(e);
        }
        Self {
            name,
            cluster,
            units,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
enum Unit {
    // The size of the gap
    G(usize),
    // Start contig,Start unit, end contig, end unit.
    E(u16, u16, u16, u16),
}

fn summarize_reads(reads: &[EncodedRead], assign: &HashMap<String, u8>) -> Vec<Read> {
    reads.iter().map(|read| Read::from(read, assign)).collect()
}
