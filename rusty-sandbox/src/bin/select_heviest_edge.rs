extern crate bio;
extern crate bio_utils;
extern crate rusty_sandbox;
use bio_utils::sam;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;


fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    // Graph construction.
    let (graph, id_to_index) = rusty_sandbox::construct_graph(&args[1])?;
    // Pairing each node with the heaviest edges.
    eprintln!("There are {} nodes.", graph.len());
    let pairs: Vec<_> = graph
        .into_iter()
        .map(|edges| {
            edges
                .into_iter()
                .fold(None, |sofar, (c_idx, c_score, c_len)| match sofar {
                    None => Some((c_idx, c_score, c_len)),
                    Some((_, score, _)) if c_score > score => Some((c_idx, c_score, c_len)),
                    Some((_, score, len)) if c_score == score && c_len > len => {
                        Some((c_idx, c_score, c_len))
                    }
                    Some(sofar) => Some(sofar),
                })
        })
        .map(|e| match e {
            Some((idx, _, _)) => Some(idx),
            None => None,
        })
        .collect();
    // Outputting.
    let (mut total_edge, mut removed_edge) = (0, 0);
    let mut wtr = BufWriter::new(std::io::stdout());
    for line in BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        if line.starts_with('@') {
            writeln!(&mut wtr, "{}", line)?;
        } else if let Some(record) = sam::Sam::new(&line) {
                total_edge += 1;
                let (q_idx, r_idx) = (id_to_index[record.q_name()], id_to_index[record.r_name()]);
                if pairs[q_idx] == Some(r_idx) || pairs[r_idx] == Some(q_idx) {
                    writeln!(&mut wtr, "{}", line)?;
                } else {
                    removed_edge += 1;
                }
            }
    }
    eprintln!("Heviest Edge selection.");
    eprintln!(
        "Total edges:{}. {} Edges were removed. {} edges remains.",
        total_edge, removed_edge, total_edge - removed_edge,
    );
    Ok(())
}
