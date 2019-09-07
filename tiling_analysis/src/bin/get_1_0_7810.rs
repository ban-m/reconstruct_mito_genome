extern crate tiling_analysis;
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
    let mut i = 0;
    for read in reads {
        let idx = read
            .iter()
            .take_while(|&unit| {
                let is_ok = match unit {
                    &Unit::Encoded(c, u, s) => c == 1 && u == 0 && (s == 7 || s == 8 || s == 10),
                    _ => false,
                };
                !is_ok
            })
            .count();
        if idx == read.len(){
            continue;
        }
        let read = EncodedRead {
            id: format!("{:2}", i),
            read: read.read.clone(),
        };
        println!("{}", read);
        i += 1;
    }
    Ok(())
}
