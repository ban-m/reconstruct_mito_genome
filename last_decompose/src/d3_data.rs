use super::find_breakpoint::{Cluster, ReadClassify};
use last_tiling::{Contigs, EncodedRead};
use serde::*;
use std::collections::HashMap;
pub fn convert_to_d3_data(
    contigs: &Contigs,
    reads: &[EncodedRead],
    clusters: &[Cluster],
) -> Summary {
    let contigs = summarize_contig(&contigs, &reads);
    let reads = summarize_reads(&reads, &clusters);
    Summary {
        contigs,
        reads,
        clusters: clusters.to_vec(),
    }
}

pub fn convert_to_d3_data_with_assign(
    contigs: &Contigs,
    reads: &[EncodedRead],
    clusters: &[Cluster],
    assignments: &HashMap<String, u8>,
) -> Summary {
    let contigs = summarize_contig(&contigs, &reads);
    let reads = summarize_reads_with_assignments(&reads, assignments);
    Summary {
        contigs,
        reads,
        clusters: clusters.to_vec(),
    }
}

pub fn convert_result_to_d3_data(
    contigs: &Contigs,
    reads: &[EncodedRead],
    assignments: &HashMap<String, u8>,
) -> ResultSummary {
    let contigs = summarize_contig(&contigs, &reads);
    let reads = summarize_reads_with_assignments(reads, assignments);
    ResultSummary { contigs, reads }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResultSummary {
    contigs: Vec<Contig>,
    reads: Vec<Read>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Summary {
    contigs: Vec<Contig>,
    reads: Vec<Read>,
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

fn summarize_contig(
    contigs: &last_tiling::Contigs,
    reads: &[last_tiling::EncodedRead],
) -> Vec<Contig> {
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
            if let last_tiling::unit::ChunkedUnit::En(encode) = unit {
                if first {
                    cs[encode.contig as usize].start_stop[encode.unit as usize] += 1;
                    first = false;
                }
                cs[encode.contig as usize].coverages[encode.unit as usize] += 1
            }
        }
        if let Some(last_tiling::unit::ChunkedUnit::En(encode)) =
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
    cluster: i32,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
enum Unit {
    // The size of the gap
    G(usize),
    E(u16, u16),
}

fn summarize_reads(reads: &[last_tiling::EncodedRead], clusters: &[Cluster]) -> Vec<Read> {
    reads
        .iter()
        .map(|read| {
            let name = read.id().to_string();
            let units = read
                .seq()
                .iter()
                .map(|e| match e {
                    last_tiling::unit::ChunkedUnit::Gap(gp) => Unit::G(gp.len()),
                    last_tiling::unit::ChunkedUnit::En(en) => Unit::E(en.contig, en.unit),
                })
                .collect();
            let cluster = get_cluster(read, clusters);
            Read {
                name,
                units,
                cluster,
            }
        })
        .collect()
}
fn summarize_reads_with_assignments(
    reads: &[EncodedRead],
    assign: &HashMap<String, u8>,
) -> Vec<Read> {
    reads
        .iter()
        .map(|read| {
            let name = read.id().to_string();
            let units = read
                .seq()
                .iter()
                .map(|e| match e {
                    last_tiling::unit::ChunkedUnit::Gap(gp) => Unit::G(gp.len()),
                    last_tiling::unit::ChunkedUnit::En(en) => Unit::E(en.contig, en.unit),
                })
                .collect();
            let cluster = assign.get(&name).map(|&e| e as i32).unwrap_or(-1);
            Read {
                name,
                units,
                cluster,
            }
        })
        .collect()
}

fn get_cluster(read: &EncodedRead, clusters: &[Cluster]) -> i32 {
    clusters
        .iter()
        .filter(|cluster| cluster.has(read.id()))
        .map(|cluster| cluster.id as i32)
        .nth(0)
        .unwrap_or(-1)
}
