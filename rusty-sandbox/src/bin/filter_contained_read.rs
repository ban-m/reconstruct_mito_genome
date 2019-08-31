extern crate bio;
extern crate bio_utils;
use bio_utils::sam;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
const THR: usize = 100;
fn id_contained_alignment(line: String) -> Option<String> {
    let sam = match sam::Sam::new(&line) {
        Some(res) => res,
        None => return None,
    };
    use bio_utils::sam::Op;
    let cigar_ref = sam.cigar();
    let is_contained = 
         match cigar_ref.first() {
            Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
            _ => 0,
        } < THR
        && match cigar_ref.last() {
            Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
            _ => 0,
        } < THR;
    if is_contained {
        Some(sam.q_name().to_string())
    } else {
        None
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    // determin contained reads.
    let contained_reads_id: HashSet<_> = {
        let input_sam = BufReader::new(File::open(&Path::new(&args[1]))?);
        input_sam
            .lines()
            .filter_map(|e| e.ok())
            .filter(|line| !line.starts_with('@'))
            .filter_map(id_contained_alignment)
            .collect()
    };
    // Filtering start.
    let mut total_edges = 0;
    let mut removed_edges = 0;
    {
        let mut output_sam = BufWriter::new(File::create(&Path::new(&args[3]))?);
        for line in BufReader::new(File::open(&Path::new(&args[1]))?)
            .lines()
            .filter_map(|e| e.ok())
        {
            if line.starts_with('@') {
                writeln!(&mut output_sam, "{}", line)?;
            } else {
                total_edges += 1;
                let sam = sam::Sam::new(&line).unwrap();
                if contained_reads_id.contains(sam.q_name()) {
                    removed_edges += 1;
                } else {
                    writeln!(&mut output_sam, "{}", line)?;
                }
            }
        }
    }
    let mut total_nodes = 0;
    let removed_nodes = contained_reads_id.len();
    {
        let mut output_fastq = bio::io::fastq::Writer::to_file(&Path::new(&args[4]))?;
        for read in bio::io::fastq::Reader::from_file(&Path::new(&args[2]))?
            .records()
            .filter_map(|e| e.ok())
        {
            total_nodes += 1;
            if !contained_reads_id.contains(read.id()) {
                output_fastq.write_record(&read)?;
            }
        }
    }
    eprintln!("Contained Read removed.");
    eprintln!(
        "Total:{} Edges. {} edges were removed. {} edges remain.",
        total_edges, removed_edges, total_edges - removed_edges
    );
    eprintln!(
        "Total:{} Nodes. {} nodes were removed. {} nodes remain.",
        total_nodes, removed_nodes,
        total_nodes - removed_nodes
    );
    Ok(())
}
