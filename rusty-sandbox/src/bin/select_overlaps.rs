extern crate bio_utils;
use bio_utils::sam;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::fs::File;
use std::path::Path;
    
const THR: usize = 100;
fn is_overlap(record: &sam::Sam, ref_length: usize) -> bool {
    use sam::Op;
    let cigar = record.cigar();
    let (ref_begin, ref_end) = {
        let (start, end) = record.get_range();
        (start, ref_length - end)
    };
    let query_begin = match cigar.first() {
        Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
        _ => 0,
    };
    let query_end = match cigar.last() {
        Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
        _ => 0,
    };
    (ref_begin < THR && query_end < THR) || (query_begin < THR && ref_end < THR)
}

fn main() -> std::io::Result<()> {
    let args:Vec<_> = std::env::args().collect();
    let mut wtr = BufWriter::new(std::io::stdout());
    let mut id_to_length = HashMap::new();
    let mut total_edges = 0;
    let mut removed_edges = 0;
    for line in BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        if line.starts_with('@') {
            if line.starts_with("@SQ") {
                let line: Vec<_> = line.split('\t').collect();
                let name = line[1].split(':').nth(1).unwrap().to_string();
                let length: usize = line[2]
                    .split(':')
                    .nth(1)
                    .and_then(|e| e.parse().ok())
                    .unwrap();
                id_to_length.insert(name, length);
            }
            writeln!(&mut wtr, "{}", line)?;
        } else if let Some((ref_length,record)) =
                sam::Sam::new(&line).map(|e| (id_to_length[e.q_name()],e))
            {
                total_edges += 1;
                if is_overlap(&record, ref_length) {
                    writeln!(&mut wtr, "{}", line)?;
                }else{
                    removed_edges += 1;
                }
            }
    }
    eprintln!("Overlap filtering done.");
    eprintln!("There are {} edges. {} Edges were removed. {} Edges remain.",
              total_edges,removed_edges, total_edges - removed_edges);
    Ok(())
}
