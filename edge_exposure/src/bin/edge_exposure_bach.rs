extern crate bio_utils;
extern crate edge_exposure;
extern crate rayon;
use edge_exposure::pick_first_cycle;
use rayon::prelude::*;
use bio_utils::sam::Sam;
use std::io::{BufRead,BufReader,BufWriter,Write};
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
    let times:usize = args[2].parse().unwrap();
    let mut log = BufWriter::new(File::create(&Path::new(&args[3]))?);
    let result:Vec<_> = (0..times).into_par_iter()
        .map(|e|(e,pick_first_cycle(nodes_size,&edges))).collect();
    for (bachnum,(picked_edges,indecies)) in result{
        writeln!(log,"{}\t{}\t{}",bachnum,picked_edges,indecies.len()).unwrap();
        for index in indecies{
            println!("{}\t{}",bachnum, ids[index]);
        }
    }
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

