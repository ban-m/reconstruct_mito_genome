extern crate rayon;
extern crate serde;
extern crate tiling_analysis;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate rand;
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use tiling_analysis::encoded_read;
use tiling_analysis::encoded_read::{EncodedRead, Unit};
use tiling_analysis::Alignment;
use std::collections::HashMap;
use rand::{thread_rng,seq::IteratorRandom};

const GAP: i64 = -4;
#[derive(Serialize, Deserialize)]
struct Graph {
    nodes: Vec<(usize, u8)>,      // Node number, its color
    edges: Vec<(usize, usize, i64)>, // (from, to, attract)
}

fn score<'r, 's>(a: &'r Unit, b: &'s Unit) -> i64 {
    match (a, b) {
        (&Unit::Encoded(u1, s1, ss1), &Unit::Encoded(u2, s2, ss2)) => {
            if u1 == u2 && s1 == s2 && ss1 == ss2{
                1
            } else {
                -2
            }
        }
        (&Unit::Gap(l1), &Unit::Gap(l2)) => {
            let diff = if l1 < l2 { l2 - l1 } else { l1 - l2 };
            let ratio = diff as f64 / l1.max(l2) as f64;
            (GAP as f64 * (1.0 - ratio)).floor() as i64
        }
        _ => -2,
    }
}

fn construct_edges<'a>(
    units: Vec<(u8, Vec<&'a EncodedRead>)>,
    thr: i64,
) -> Vec<(&'a str, &'a str, i64)> {
    units
        .into_par_iter()
        .flat_map(|(_, reads_contains_unit)| {
            let len = reads_contains_unit.len();
            let mut res = vec![];
            for i in 0..len {
                for j in i + 1..len {
                    let a = reads_contains_unit[i];
                    let b = reads_contains_unit[j];
                    let alignment = Alignment::new(a, b, GAP, score);
                    if alignment.score >= thr {
                        res.push((a.id(), b.id(), alignment.score));
                    }
                }
            }
            res
        })
        .collect()
}

fn construct_graph_inner(
    records: Vec<EncodedRead>,
    names: HashMap<String, usize>,
    thr: i64,
) -> Graph {
    let units = encoded_read::collect_units(&records);
    let nodes: Vec<_> = records
        .iter()
        .filter_map(|read| Some((names[read.id()], read.color()?)))
        .collect();
    let edges = construct_edges(units, thr)
        .into_iter()
        .map(|(a, b, score)| (names[a], names[b], score))
        .collect();
    Graph { nodes, edges }
}

fn construct_graph(file: &str, thr: i64) -> std::io::Result<Graph> {
    let records: Vec<_> = BufReader::new(std::fs::File::open(&std::path::Path::new(file))?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| EncodedRead::new(&e))
        .collect();
    // To conpute the full graph, remove two line below.
    let mut rng = thread_rng();
    let records:Vec<_> = records.into_iter().choose_multiple(&mut rng, 500);
    let names: HashMap<String, usize> = records
        .iter()
        .map(|e| e.id().to_string())
        .enumerate()
        .map(|(idx, id)| (id, idx))
        .collect();
    Ok(construct_graph_inner(records, names, thr))
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let thr: i64 = args[2].parse().unwrap();
    let result = construct_graph(&args[1], thr)?;
    println!("{}", serde_json::ser::to_string(&result).unwrap());
    Ok(())
}
