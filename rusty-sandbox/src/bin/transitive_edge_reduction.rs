extern crate bio;
extern crate bio_utils;
extern crate rusty_sandbox;
use bio_utils::sam;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

type Graph = Vec<HashSet<usize>>;
fn construct_graph(path: &str) -> std::io::Result<(Graph, HashMap<String, usize>)> {
    let (graph, id_to_index) = rusty_sandbox::construct_graph(path)?;
    let graph:Vec<_> = graph.into_iter().map(|vs| vs.into_iter().map(|e|e.0).collect()).collect();
    Ok((graph, id_to_index))
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    // Graph construction.
    let (mut graph, id_to_index) = construct_graph(&args[1])?;
    // Transitive edge reduction.
    for idx in 0..graph.len() {
        // for each read, check whether the grandchildren would overlap with child.
        let mut removed_child = HashSet::new();
        for child in &graph[idx] {
            if !removed_child.contains(child) {
                let grand_children = &graph[*child];
                for shared_child in grand_children.intersection(&graph[idx]) {
                    removed_child.insert(*shared_child);
                }
            }
        }
        graph[idx] = graph[idx].difference(&removed_child).copied().collect();
    }
    // Output the graph.
    let (mut total_edge, mut removed_edge) = (0, 0);
    let mut wtr = BufWriter::new(std::io::stdout());
    for line in BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        if line.starts_with('@') {
            writeln!(&mut wtr, "{}", line)?;
        } else  if let Some(record) = sam::Sam::new(&line) {
                total_edge += 1;
                let (q_idx, r_idx) = (id_to_index[record.q_name()], id_to_index[record.r_name()]);
                if graph[q_idx].contains(&r_idx) || graph[r_idx].contains(&q_idx) {
                    writeln!(&mut wtr, "{}", line)?;
                } else {
                    removed_edge += 1;
                }
        }
    }
    eprintln!(
        "Total edges:{}. {} Edges were removed. {} Edges remain.",
        total_edge, removed_edge, total_edge - removed_edge
    );
    Ok(())
}
