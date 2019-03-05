extern crate bio_utils;
extern crate edge_exposure;
use bio_utils::sam::Sam;
use std::io::{BufRead,BufReader};
use std::fs::File;
use std::path::Path;
fn main() -> std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let (sam_records,ids):(Vec<_>,Vec<_>) = BufReader::new(
        &File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e|e.ok())
        .partition(|e|e.starts_with('@'));
    let mut ids:Vec<String> = ids.into_iter().filter(|e|e.starts_with("@SN"))
        .filter_map(|e|e.split('\t').nth(1)
                    .map(|e|e.trim_start_matches("SN:").to_string()))
        .collect();
    let sam_records:Vec<_> = sam_records.into_iter()
        .filter_map(|e|Sam::new(&e))
        .collect();
    ids.sort();
    let edges = construct_edges(&sam_records,&ids);
    let nodes_size = ids.len();
    eprintln!("There are {} nodes and {} edges",nodes_size, edges.len());
    
    Ok(())
}
fn construct_edges(sam_records:&[Sam],ids:&[String])->Vec<(usize,usize)>{
    sam_records.iter()
        .map(|e|(to_idx(e.q_name(),ids),
                 to_idx(e.r_name(),ids)))
        .collect()
}

fn to_idx(name:&str,ids:&[String])->usize{
    let name = name.to_string();
    ids.binary_search_by(|probe|probe.cmp(&name)).unwrap()
}
