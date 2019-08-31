extern crate bio;
extern crate bio_utils;
extern crate rand;
extern crate rayon;
extern crate rusty_sandbox;
use rayon::prelude::*;
use std::collections::{HashMap};
use std::io::{BufWriter, Write};
use std::path::Path;

#[inline]
fn is_different_base(x: u8, y: u8) -> bool {
    match (x, y) {
        (b'a', b'a')
        | (b'A', b'a')
        | (b'a', b'A')
        | (b'A', b'A')
        | (b'c', b'c')
        | (b'C', b'c')
        | (b'c', b'C')
        | (b'C', b'C')
        | (b'g', b'g')
        | (b'G', b'g')
        | (b'g', b'G')
        | (b'G', b'G')
        | (b't', b't')
        | (b'T', b't')
        | (b't', b'T')
        | (b'T', b'T') => false,
        _ => true,
    }
}

fn error_profile(
    alignments: Vec<bio_utils::sam::Sam>,
    queries: &HashMap<String, bio::io::fastq::Record>,
    reference: &bio::io::fastq::Record,
) -> Vec<(usize, usize, usize, usize)> {
    let reference = reference.seq();
    let mut error_profile: Vec<(usize, usize, usize, usize)> = vec![(0, 0, 0, 0); reference.len()];
    for record in alignments {
        let query = if record.is_template() {
            queries[record.q_name()].seq().to_vec()
        } else {
            bio::alphabets::dna::revcomp(queries[record.q_name()].seq())
        };
        let mut ref_pos = record.pos() - 1;
        let mut query_pos = 0;
        use bio_utils::sam::Op;
        for op in record.cigar() {
            match op {
                Op::Align(l) | Op::Match(l) | Op::Mismatch(l) => {
                    for i in 0..l {
                        if is_different_base(query[query_pos + i], reference[ref_pos + i]) {
                            error_profile[ref_pos + i].0 += 1;
                        }
                        error_profile[ref_pos + i].3 += 1;
                    }
                    query_pos += l;
                    ref_pos += l;
                }
                Op::Deletion(l) => {
                    for i in 0..l {
                        error_profile[ref_pos + i].1 += 1;
                    }
                    ref_pos += l;
                }
                Op::Insertion(l) => {
                    error_profile[ref_pos].2 += l;
                    query_pos += l;
                }
                Op::SoftClip(l) | Op::HardClip(l) => query_pos += l,
                _ => {}
            };
        }
    }
    error_profile
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("[sam file] [query file] [out dir]");
        return Ok(());
    }
    let sam_records = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    let queries: HashMap<_,_> = bio::io::fastq::Reader::from_file(&Path::new(&args[2]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e|(e.id().to_string(),e))
        .collect();
    // Error Profiles
    sam_records
        .into_par_iter()
        .filter_map(|(name,alignments)| {
            eprintln!("Profilling {}", name);
            let reference = &queries[&name];
            let error_profile = error_profile(alignments, &queries, &reference);
            let mut wtr = BufWriter::new(
                std::fs::File::create(&Path::new(&format!(
                    "{}/{}",
                    &args[3],
                    name.replace('/', "")
                )))
                .ok()?,
            );
            writeln!(&mut wtr, "Index\tSubst\tDel\tIns\tDepth").ok()?;
            for (idx, (subst, del, ins, depth)) in error_profile.into_iter().enumerate() {
                writeln!(&mut wtr, "{}\t{}\t{}\t{}\t{}", idx, subst, del, ins, depth).ok()?;
            }
            Some(())
        })
        .count();
    Ok(())
}
