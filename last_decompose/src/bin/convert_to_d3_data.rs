extern crate last_tiling;
#[macro_use]
extern crate serde;
extern crate last_decompose;
extern crate serde_json;
use last_decompose::find_breakpoint::ReadClassify;
use std::io::{BufReader, BufWriter};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let contigs: last_tiling::Contigs =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let reads: Vec<last_tiling::EncodedRead> =
        serde_json::de::from_reader(std::fs::File::open(&args[2]).map(BufReader::new)?).unwrap();
    let repeats: Vec<last_tiling::repeat::RepeatPairs> = last_tiling::repeat::open(&args[3])?;
    let cr = {
        let reads: Vec<_> = reads
            .iter()
            .map(|read| last_decompose::ERead::new(read.clone()))
            .collect();
        last_decompose::critical_regions(&reads, &contigs, &repeats)
    };
    let contigs = summarize_contig(&contigs, &reads);
    let reads = summarize_reads(&reads, &cr);
    let summary = Summary {
        contigs,
        reads,
        critical_regions: cr,
    };
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    serde_json::ser::to_writer(&mut stdout, &summary).unwrap();
    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Summary {
    contigs: Vec<Contig>,
    reads: Vec<Read>,
    critical_regions: Vec<last_decompose::CriticalRegion>,
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
            match unit {
                last_tiling::unit::ChunkedUnit::En(encode) => {
                    if first {
                        cs[encode.contig as usize].start_stop[encode.unit as usize] += 1;
                        first = false;
                    }
                    cs[encode.contig as usize].coverages[encode.unit as usize] += 1
                }
                _ => {}
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
    Gap(usize),
    Encode(u16, u16),
}

fn summarize_reads(
    reads: &[last_tiling::EncodedRead],
    cr: &[last_decompose::CriticalRegion],
) -> Vec<Read> {
    reads
        .iter()
        .map(|read| Read {
            name: read.id().to_string(),
            units: read
                .seq()
                .into_iter()
                .map(|e| match e {
                    last_tiling::unit::ChunkedUnit::Gap(gp) => Unit::Gap(gp.len()),
                    last_tiling::unit::ChunkedUnit::En(en) => Unit::Encode(en.contig, en.unit),
                })
                .collect(),
            cluster: get_cluster(read, cr),
        })
        .collect()
}

fn get_cluster(read: &last_tiling::EncodedRead, cr: &[last_decompose::CriticalRegion]) -> i32 {
    let read = last_decompose::ERead::new(read.clone());
    cr.iter()
        .enumerate()
        .filter_map(|(idx, cr)| {
            if cr.along_with(&read) {
                Some(idx)
            } else {
                None
            }
        })
        .nth(0)
        .map(|idx| idx as i32)
        .unwrap_or(-1)
}
