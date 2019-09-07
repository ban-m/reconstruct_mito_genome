extern crate tiling_analysis;
use std::collections::{HashMap, HashSet};
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
    let nodes_kind: HashSet<_> = reads.iter().flat_map(|e| e.iter()).collect();
    let total_edges: usize = reads.iter().map(|e| e.len() - 1).sum();
    eprintln!(
        "There are {} edges. {} nodes.",
        total_edges,
        total_edges + reads.len()
    );
    eprintln!("There are {} reads", reads.len());
    eprintln!("There are {} kinds of nodes", nodes_kind.len());
    eprintln!("Enumerating non canonical edges.");

    let mut non_canonical_edges: HashMap<_, usize> = HashMap::new();
    for read in reads {
        for seq in read.read.windows(2) {
            if !is_canonical(&seq[0], &seq[1]) {
                let entry = non_canonical_edges
                    .entry((seq[0].clone(), seq[1].clone()))
                    .or_default();
                *entry += 1;
            }
        }
    }
    for ((u, v), count) in non_canonical_edges {
        if count >= 3 {
            println!("{},{},{}", u, v, count);
        }
    }
    Ok(())
}

fn is_canonical(u: &Unit, v: &Unit) -> bool {
    match (u, v) {
        (&Unit::Encoded(c1, u1, s1), &Unit::Encoded(c2, u2, s2)) => {
            (c1 == c2 && u1 == u2 && s1 + 1 == s2)
                || (c1 == c2 && u1 == u2 && s1 + 2 == s2)
                || (c1 == c2 && u1 != u2 && s2 == 0)
                || (c1 == c2 && u1 == u2 && s1 - 1 == s2)
                || (c1 == c2 && u1 == u2 && s1 - 2 == s2)
                || (c1 == c2 && u1 != u2 && s2 == 0)
        }
        (&Unit::Gap(_), _) => true,
        (_, &Unit::Gap(_)) => true,
    }
}
