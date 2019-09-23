extern crate tiling_analysis;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tiling_analysis::encoded_read::*;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| EncodedRead::new(&e))
        .filter(|e| !e.is_empty())
        .collect();
    let unit_count:HashMap<_,usize> =
        reads
            .into_iter()
            .flat_map(|e| e.read)
            .fold(HashMap::new(), |mut res, unit| {
                let entry = res.entry(unit).or_default();
                *entry += 1;
                res
            });
    println!("ctg\tunit\tsubunit\tcount");
    for (unit, count) in unit_count {
        if let Unit::Encoded(ctg, unit, subunit) = unit {
            println!("{}\t{}\t{}\t{}", ctg, unit, subunit, count);
        }
    }
    Ok(())
}
