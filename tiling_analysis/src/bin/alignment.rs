extern crate rayon;
extern crate tiling_analysis;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tiling_analysis::encoded_read::*;
use tiling_analysis::Alignment;
const GAP:i64  = -3;
fn score<'r, 's>(a: &'r Unit, b: &'s Unit) -> i64 {
    match (a, b) {
        (&Unit::Encoded(u1, s1), &Unit::Encoded(u2, s2)) => {
            if u1 == u2 && s1 == s2 {
                20
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

use std::collections::HashSet;
fn collect_units<'a>(reads: &'a Vec<EncodedRead>) -> Vec<(usize, Vec<&'a EncodedRead>)> {
    let units: HashSet<_> = reads
        .iter()
        .flat_map(|e| e.iter().filter_map(|e| e.get_unit()))
        .collect();
    units
        .iter()
        .map(|&unit| {
            (
                unit,
                reads.iter().filter(|read| read.contains(unit)).collect(),
            )
        })
        .collect()
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    eprintln!("Opening and parsing...");
    let thr:i64 = args[2].parse().unwrap();
    let reads: Vec<_> = BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| EncodedRead::new(&e))
        .filter(|e| !e.is_empty())
        .collect();
    eprintln!("There are {} reads", reads.len());
    let units = collect_units(&reads);
    let result:Vec<_> = units
        .into_par_iter()
        .map(|(unit, reads_contains_unit)| {
            let start = std::time::Instant::now();
            let len = reads_contains_unit.len();
            let mut res = String::new();
            for i in 0..len {
                for j in i+1..len {
                    let alignment =
                        Alignment::new(reads_contains_unit[i], reads_contains_unit[j], GAP, score);
                    if alignment.score >= thr {
                        res.push_str(&alignment.to_string());
                    }
                }
            }
            let end = std::time::Instant::now();
            eprintln!("Unit {} contains {}.... Aligned in {:?}.", unit, len, end - start);
            res
        })
        .collect();
    for r in result {
        println!("{}", r);
    }
    Ok(())
}
