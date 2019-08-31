extern crate bio;
extern crate bio_utils;
extern crate rust_htslib;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
const THR: usize = 100;
use bio_utils::sam;
pub fn reverse_index(reader: &rust_htslib::bam::Reader) -> Vec<String> {
    let header = reader.header();
    let mut res = vec![];
    for name in header.target_names() {
        let idx = header.tid(name).unwrap() as usize;
        (res.len()..=idx).for_each(|_| res.push(String::new()));
        res[idx] = String::from_utf8_lossy(name).to_string();
    }
    res
}

pub fn open_sam_into_hashmap(
    path: &str,
) -> std::io::Result<HashMap<String, Vec<bio_utils::sam::Sam>>> {
    let mut result = HashMap::new();
    let mut buf = fs::File::open(&std::path::Path::new(path)).map(BufReader::new)?;
    let mut line = String::new();
    while buf.read_line(&mut line)? > 0 {
        if let Some(sam) = bio_utils::sam::Sam::new(&line) {
            let entry = result
                .entry(sam.r_name().to_string())
                .or_insert_with(Vec::new);
            entry.push(sam);
        }
        line.clear();
    }
    Ok(result)
}

pub fn open_sam_into_vec(path: &str) -> std::io::Result<Vec<bio_utils::sam::Sam>> {
    let mut result = vec![];
    let mut buf = fs::File::open(&std::path::Path::new(path)).map(BufReader::new)?;
    let mut line = String::new();
    while buf.read_line(&mut line)? > 0 {
        if let Some(sam) = bio_utils::sam::Sam::new(&line) {
            result.push(sam);
        }
    }
    Ok(result)
}

pub fn open_fastq_into_hashmap(
    path: &str,
) -> std::io::Result<HashMap<String, bio::io::fastq::Record>> {
    bio::io::fastq::Reader::from_file(&std::path::Path::new(&path)).map(|rec| {
        rec.records()
            .filter_map(|e| e.ok())
            .map(|e| (e.id().to_string(), e))
            .collect()
    })
}

/// Parse the header "@SQ" field in sam file.
pub fn extract_id_and_length(line: &str) -> (String, usize) {
    let line: Vec<_> = line.split('\t').collect();
    let name = line[1].split(':').nth(1).unwrap().to_string();
    let length: usize = line[2]
        .split(':')
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap();
    (name, length)
}

pub fn determine_direction(record: &sam::Sam, ref_length: usize) -> Option<(&str, &str)> {
    // return (from to).
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
    if ref_begin < THR && query_end < THR {
        // query to reference.
        Some((record.q_name(), record.r_name()))
    } else if query_begin < THR && ref_end < THR {
        // ref to query.
        Some((record.r_name(), record.q_name()))
    } else {
        // panic!("Something wrong. You might've missed the overlap detection step.");
        None
    }
}

type Graph = Vec<Vec<(usize, usize, usize)>>;
/// Construct a directed graph from specified SAM file.
/// Each node list is (to, mapping quality, alignment length).
pub fn construct_graph(path: &str) -> std::io::Result<(Graph, HashMap<String, usize>)> {
    let mut id_to_index = HashMap::new();
    let mut id_to_length = HashMap::new();
    let mut node_number = 0;
    // graph[i] = [(j,q,l)] <=> there is edges from i to j, with quality score q, and alignment length l.
    let mut graph = vec![];
    for line in BufReader::new(File::open(&Path::new(path))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        if line.starts_with("@SQ") {
            // new node.
            let (id, length) = extract_id_and_length(&line);
            id_to_index.insert(id.clone(), node_number);
            id_to_length.insert(id, length);
            node_number += 1;
        } else if !line.starts_with('@') {
            let sam = match sam::Sam::new(&line) {
                Some(res) => res,
                None => continue,
            };
            let length = match id_to_length.get(sam.q_name()) {
                Some(len) => len,
                None => continue,
            };
            let (from, to) = match determine_direction(&sam, *length) {
                Some(res) => res,
                None => continue,
            };
            let from = match id_to_index.get(from) {
                Some(res) => *res,
                None => continue,
            };
            let to = match id_to_index.get(to) {
                Some(res) => *res,
                None => continue,
            };
            (graph.len()..=from).for_each(|_| graph.push(vec![]));
            let (start, end) = sam.mapped_region();
            graph[from].push((to, sam.mapq(), end - start));
        }
    }
    (graph.len()..node_number).for_each(|_| graph.push(vec![]));
    Ok((graph, id_to_index))
}

/// Take the adjacency list, return the connected component.
pub fn calc_connected_component(graph: &[Vec<usize>]) -> Vec<HashSet<usize>> {
    let mut count = 0;
    let mut has_assigned = vec![false; graph.len()];
    let mut preorder = vec![None; graph.len()];
    let mut edge_index = vec![0; graph.len()];
    let mut idx = 0;
    let mut result = vec![];
    while idx < graph.len() {
        let mut dfs = vec![];
        dfs.push(idx);
        let mut not_assined_nodes = vec![];
        let mut intermidiate_nodes = vec![];
        'outer: while !dfs.is_empty() {
            let v = *dfs.last().unwrap(); // safe.
            if preorder[v].is_none() {
                // This is preorder execution.
                not_assined_nodes.push(v);
                intermidiate_nodes.push(v);
                preorder[v] = Some(count);
                count += 1;
            }
            while edge_index[v] < graph[v].len() {
                // There is still an edge to be used.
                let u = graph[v][edge_index[v]];
                edge_index[v] += 1;
                if preorder[u].is_none() {
                    // And that node has never been arrived.
                    dfs.push(u);
                    continue 'outer; // Early drop, continue to next loop.
                } else if !has_assigned[u] {
                    let order = preorder[u].unwrap(); // safe.
                    while let Some(top) = intermidiate_nodes.pop() {
                        if top <= order {
                            intermidiate_nodes.push(top);
                            break;
                        }
                    }
                }
            }
            // Here, no nodes was added in the while loop above.
            // Continue to postoder procedure.
            let v = dfs.pop().unwrap();
            let top = intermidiate_nodes.pop().unwrap();
            if top == v {
                let mut components = HashSet::new();
                while let Some(res) = not_assined_nodes.pop() {
                    components.insert(res);
                    has_assigned[res] = true;
                    if res == v {
                        break;
                    }
                }
                result.push(components);
            } else {
                intermidiate_nodes.push(top);
            }
        }
        // Seach for next tree.
        while idx < graph.len() && preorder[idx].is_some() {
            idx += 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn work() {}
    #[test]
    fn directed_component() {
        let graph = vec![
            vec![1],
            vec![2],
            vec![3, 4],
            vec![0],
            vec![5],
            vec![6],
            vec![4],
        ];
        let dc1: HashSet<usize> = vec![0, 1, 2, 3].into_iter().collect();
        let dc2: HashSet<usize> = vec![4, 5, 6].into_iter().collect();
        let result = calc_connected_component(&graph);
        eprintln!("{:?}", result);
        assert_eq!(result.len(), 2);
        assert!(result[0] == dc1 || result[0] == dc2);
        assert!(result[1] == dc1 || result[1] == dc2);
    }
}
