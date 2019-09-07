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
    eprintln!("There are {} reads", reads.len());
    for read in reads {
        for w in read.read.windows(2) {
            match w {
                [Unit::Encoded(c1, u1, s1), Unit::Encoded(c2, u2, s2)] => {
                    let good = (c1 == c2 && u1 == u2 && s1 + 1 == *s2)
                        || (c1 == c2 && u1 != u2 && *s2 == 0)
                        || (c1 == c2 && u1 == u2 && s1 - 1 == *s2)
                        || (c1 == c2 && u1 != u2 && *s2 == 0);
                    if !good {
                        println!("{}", read.to_output());
                        break;
                    }
                }
                _ => continue,
            }
        }
    }
    Ok(())
}
