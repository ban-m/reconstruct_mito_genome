extern crate bio;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
#[derive(Serialize, Deserialize, Debug, Clone)]
struct RepeatAnnotation {
    repeats: Vec<Repeat>,
}

// impl RepeatAnnotation {
//     fn to_tidy(mut self) -> Self {
//         if self.repeats.is_empty() {
//             self
//         } else {
//             self.repeats = Self::to_tidy_inner(self.repeats);
//             self
//         }
//     }
//     // To merge consective kmer match
//     fn to_tidy_inner(mut repeats: Vec<Repeat>) -> Vec<Repeat> {
//         repeats.sort();
//         for repeat in &repeats{
//             eprintln!("{:?}", repeat);
//         }
//         let mut current = repeats[0].clone();
//         let mut predicted_next = current.predict();
//         let mut result = vec![];
//         for repeat in repeats.iter().skip(1) {
//             if &predicted_next == repeat {
//                 current.k += 1;
//             } else {
//                 result.push(current);
//                 current = repeat.clone();
//             }
//             predicted_next = repeat.predict();
//         }
//         result.push(current);
//         println!("");
//         result
//     }
// }

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Repeat {
    pos1: Position,
    pos2: Position,
    k: usize,
}

impl Repeat {
    // fn predict(&self) -> Repeat {
    //     assert!(self.pos1.is_forward);
    //     let start1 = self.pos1.start + 1;
    //     let start2 = if self.pos2.is_forward {
    //         self.pos2.start + 1
    //     } else {
    //         self.pos2.start - 1
    //     };
    //     Repeat {
    //         pos1: Position {
    //             start: start1,
    //             is_forward: self.pos1.is_forward,
    //         },
    //         pos2: Position {
    //             start: start2,
    //             is_forward: self.pos2.is_forward,
    //         },
    //         k: self.k,
    //     }
    // }

    fn new(c1: &str, pos1: usize, isf1: bool, c2: &str, pos2: usize, isf2: bool, k: usize) -> Self {
        Repeat {
            pos1: Position {
                contig: c1.to_string(),
                start: pos1,
                is_forward: isf1,
            },
            pos2: Position {
                contig: c2.to_string(),
                start: pos2,
                is_forward: isf2,
            },
            k: k,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Position {
    contig: String,
    start: usize,
    is_forward: bool,
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let repeats = parse_mummer(&args[1])?;
    let repeats = serde_json::ser::to_string_pretty(&repeats).unwrap();
    println!("{}", repeats);
    Ok(())
}

fn parse_mummer(file: &str) -> std::io::Result<RepeatAnnotation> {
    let mut buffer = vec![]; // bucket per contig
    let mut repeats = vec![]; //result
    for line in BufReader::new(File::open(&Path::new(file))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        if line.starts_with('>') {
            repeats.extend(parse_mummer_bucket(&buffer));
            buffer.clear();
        }
        buffer.push(line);
    }
    repeats.extend(parse_mummer_bucket(&buffer));
    Ok(RepeatAnnotation { repeats })
}

fn parse_mummer_bucket(buffer: &Vec<String>) -> Vec<Repeat> {
    if buffer.is_empty() {
        return vec![];
    };
    let contig1 = parse_contig_name(&buffer[0]);
    let is_reverse = !buffer[0].contains("Reverse");
    buffer
        .iter()
        .skip(1)
        .filter_map(|e| to_repeat(&contig1, is_reverse, e))
        .collect()
}

fn parse_contig_name(line: &str) -> String {
    line.trim_start_matches('>')
        .trim_end_matches("Reverse")
        .trim()
        .to_string()
}

fn to_repeat(refcontig: &str, is_forward: bool, line: &str) -> Option<Repeat> {
    let line_num = line.split_whitespace().count();
    if line_num == 3 {
        to_repeat_3fileds(refcontig, is_forward, line)
    } else if line_num == 4 {
        to_repeat_4fileds(refcontig, is_forward, line)
    } else {
        unreachable!();
    }
}

fn to_repeat_4fileds(refcontig: &str, is_forward: bool, line: &str) -> Option<Repeat> {
    let line: Vec<_> = line.split_whitespace().collect();
    let query = line[0].to_string();
    let line: Vec<usize> = line
        .into_iter()
        .skip(1)
        .filter_map(|e| e.parse().ok())
        .collect();
    if line.len() != 3 || (line[0] == line[1] && line[0] == 1) {
        None
    } else {
        // To convert 0-based
        Some(Repeat::new(
            refcontig,
            line[0] - 1,
            is_forward,
            &query,
            line[1] - 1,
            true,
            line[2],
        ))
    }
}

fn to_repeat_3fileds(refcontig: &str, is_forward: bool, line: &str) -> Option<Repeat> {
    // The repeats resides in the same contig.
    let line: Vec<usize> = line
        .split_whitespace()
        .filter_map(|e| e.parse().ok())
        .collect();
    if line.len() != 3 || (line[0] == line[1] && line[0] == 1) {
        None
    } else {
        // To convert 0-based
        Some(Repeat::new(
            refcontig,
            line[0] - 1,
            is_forward,
            refcontig,
            line[1] - 1,
            true,
            line[2],
        ))
    }
}

// fn find_repeats(path: &str, k: usize) -> std::io::Result<Vec<RepeatAnnotation>> {
//     let input = bio::io::fasta::Reader::from_file(&std::path::Path::new(path))?;
//     Ok(input
//         .records()
//         .filter_map(|e| e.ok())
//         .map(|contig| find_repeats_from_contig(contig, k))
//         .collect())
// }

// use std::collections::HashMap;
// fn find_repeats_from_contig(contig: bio::io::fasta::Record, k: usize) -> RepeatAnnotation {
//     let kmers = collect_kmers(contig.seq(), k);
//     let name = contig.id().to_string();
//     let repeats: Vec<_> = contig
//         .seq()
//         .windows(k)
//         .enumerate()
//         .flat_map(|(p, kmer)| find_repeats_by_kmer(p, &kmers[kmer.as_ref()], k))
//         .collect();

//     RepeatAnnotation {
//         contig: name,
//         repeats: repeats,
//     }.to_tidy()
// }

// fn find_repeats_by_kmer(position: usize, kmers: &Vec<(usize, bool)>, k: usize) -> Vec<Repeat> {
//     if kmers.len() == 1 {
//         vec![]
//     } else {
//         kmers
//             .iter()
//             .filter(|(p, _)| p > &position)
//             .map(|(p, d)| Repeat::new(position, true, *p, *d, k))
//             .collect()
//     }
// }

// fn collect_kmers(seq: &[u8], k: usize) -> HashMap<Vec<u8>, Vec<(usize, bool)>> {
//     let mut result = HashMap::new();
//     let len = seq.len();
//     for (kmer, index, dir) in seq
//         .windows(k)
//         .enumerate()
//         .map(|(index, kmer)| (kmer.to_vec(), index, true))
//         .chain(
//             bio::alphabets::dna::revcomp(seq)
//                 .windows(k)
//                 .enumerate()
//                 .map(|(index, kmer)| (len - index - k, kmer))
//                 .map(|(index, kmer)| (kmer.to_vec(), index, false)),
//         )
//     {
//         let entry = result.entry(kmer).or_insert(vec![]);
//         entry.push((index, dir));
//     }
//     result
// }

// #[cfg(test)]
// pub mod tests {
//     use super::*;
//     #[test]
//     fn to_tidy_test() {
//         let repeats = vec![
//             Repeat::new(10, true, 12, true, 8),
//             Repeat::new(11, true, 13, true, 8),
//             Repeat::new(12, true, 14, true, 8),
//         ];
//         assert_eq!(
//             RepeatAnnotation::to_tidy_inner(repeats),
//             vec![Repeat::new(10, true, 12, true, 10)]
//         );
//         let repeats = vec![
//             Repeat::new(12, true, 14, true, 8),
//             Repeat::new(11, true, 13, true, 8),
//             Repeat::new(10, true, 12, true, 8),
//         ];
//         assert_eq!(
//             RepeatAnnotation::to_tidy_inner(repeats),
//             vec![Repeat::new(10, true, 12, true, 10)]
//         );
//         let repeats = vec![
//             Repeat::new(13, true, 14, true, 8),
//             Repeat::new(11, true, 13, true, 8),
//             Repeat::new(10, true, 12, true, 8),
//         ];
//         assert_eq!(
//             RepeatAnnotation::to_tidy_inner(repeats),
//             vec![Repeat::new(10, true, 12, true, 9),
//                  Repeat::new(13, true, 14, true, 8)
//             ]
//         );
//         let repeats = vec![
//             Repeat::new(12, true, 15, true, 12),
//             Repeat::new(12, true, 21, true, 12),
//             Repeat::new(13, true, 16, true, 12),
//             Repeat::new(13, true, 22, true, 12),
//             Repeat::new(12, true, 312, true, 12),
//         ];
//         assert_eq!(
//             RepeatAnnotation::to_tidy_inner(repeats),
//             vec![Repeat::new(12, true, 16, true, 13),
//                  Repeat::new(12, true, 22, true, 13),
//                  Repeat::new(12, true, 312, true, 12)
//             ]
//         );
//         let repeats = vec![
//             Repeat::new(18, true, 2, false, 11),
//             Repeat::new(19, true, 1, false, 11),
//             Repeat::new(20, true, 0, false, 11)
//         ];
//         assert_eq!(RepeatAnnotation::to_tidy_inner(repeats),
//                    vec![Repeat::new(18, true, 0,false, 13)]);

//         let repeats = vec![
//             Repeat::new(19, true, 1, false, 11),
//             Repeat::new(18, true, 2, false, 11),
//             Repeat::new(20, true, 0, false, 11)
//         ];
//         assert_eq!(RepeatAnnotation::to_tidy_inner(repeats),
//                    vec![Repeat::new(18, true, 0, false,13)]);

//     }
//     #[test]
//     fn to_predict_next() {
//         let before = Repeat::new(10, true, 20, true, 8);
//         assert_eq!(
//             Repeat::new(11, true, 21, true, 8),
//             before.predict());
//         let before = Repeat::new(10,true,19,false,9);
//         assert_eq!(
//             Repeat::new(11,true,18,false,9),
//             before.predict());
//     }

// }
