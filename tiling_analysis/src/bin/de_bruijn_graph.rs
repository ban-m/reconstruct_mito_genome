extern crate rayon;
extern crate tiling_analysis;
// use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use tiling_analysis::encoded_read::*;
const K: usize = 3;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| EncodedRead::new(&e))
        .filter(|e| !e.is_empty())
        .map(canonicalize)
        .map(gap_filling)
        .collect();
    let kmer_count: HashMap<_, u8> = reads.iter().fold(HashMap::new(), |mut res, read| {
        for kmer in read.read.windows(K) {
            let entry = res.entry(kmer.to_vec()).or_default();
            *entry += 1;
        }
        res
    });
    let filtered_kmer: HashSet<_> = kmer_count
        .into_iter()
        .filter_map(|(key, val)| if val > 2 { Some(key) } else { None })
        .collect();
    let mut edges: HashMap<_, Vec<_>> = HashMap::new();
    for read in reads.into_iter().filter(|read| read.len() > K) {
        let mut from = read.read[..K].to_vec();
        for unit in read.iter().skip(K) {
            let next = {
                let mut x = from[1..].to_vec();
                x.push(unit.clone());
                x
            };
            if filtered_kmer.contains(&from) && filtered_kmer.contains(&next) {
                let entry = edges.entry(from).or_default();
                entry.push(next.clone());
            }
            from = next;
        }
    }
    {
        // Distribution of degree
        let mut wtr = BufWriter::new(File::create("dgree.txt")?);
        for (_, v) in edges.iter() {
            writeln!(&mut wtr, "{}", v.len())?;
        }
    }
    let edge_size: usize = edges.iter().map(|(_, v)| v.len()).sum();
    eprintln!(
        "there are {} nodes, {} edges. {} are non-canonical.",
        filtered_kmer.len(),
        edge_size,
        edge_size - filtered_kmer.len()
    );
    let mut wtr = BufWriter::new(File::create("chunked_degree.txt")?);
    for (kmer, mut connecteds) in edges {
        connecteds.sort();
        connecteds.dedup();
        writeln!(&mut wtr, "{}", connecteds.len())?;
        if connecteds.len() > 1 {
            print(&kmer);
            print!("===>");
            for node in connecteds {
                print(&node);
                print!("\t");
            }
            println!()
        }
    }
    Ok(())
}

fn print(kmer: &[Unit]) {
    assert!(kmer.len() >= 1);
    for k in kmer {
        print!("{} ", k);
    }
}

fn canonicalize(read: EncodedRead) -> EncodedRead {
    let (forward, backward) =
        read.read
            .windows(2)
            .fold((0, 0), |(mut forward, mut backward), window| {
                if let [Unit::Encoded(c1, u1, s1), Unit::Encoded(c2, u2, s2)] = window {
                    if c1 == c2 {
                        // forward
                        if (u1 == u2 && s1 + 1 == *s2) || (u1 + 1 == *u2 && *s2 == 0) {
                            forward += 1;
                        }
                        if (u1 == u2 && *s1 == *s2 + 1) || (*u1 == *u2 + 1 && *s1 == 0) {
                            backward += 1;
                        }
                    }
                }
                (forward, backward)
            });
    if forward > backward {
        // eprintln!("{} is forward", read);
        read
    } else {
        // eprintln!("{} is backward", read);
        EncodedRead {
            id: read.id.to_string(),
            read: read.read.into_iter().rev().collect(),
        }
    }
}

fn gap_filling(mut read: EncodedRead) -> EncodedRead {
    // eprint!("{}->", read);
    for i in 1..read.len() - 1 {
        if let (&Unit::Encoded(c1, u1, s1), &Unit::Gap(_), &Unit::Encoded(c2, u2, s2)) =
            (&read[i - 1], &read[i], &read[i + 1])
        {
            if c1 == c2 && u1 == u2 && s1 + 2 == s2 {
                read.read[i] = Unit::Encoded(c1, u1, s1 + 1);
            } else {
                read.read[i] = Unit::Gap(100);
            }
        }
    }
    // eprintln!("{}", read);
    read
}
