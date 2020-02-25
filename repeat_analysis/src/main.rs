extern crate bio_utils;
extern crate last_tiling;
extern crate rayon;
use bio_utils::fasta;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
const IS_CURCULAR: bool = true;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reference = fasta::parse_into_vec(&args[1])?;
    let alignment: Vec<_> = last_tiling::parse_tab_file(&args[2])?;
    let parse_line = |line: String| {
        let line: Vec<&str> = line.split('\t').collect();
        let id = line[0].to_string();
        let locs: Vec<usize> = line[1..].iter().filter_map(|e| e.parse().ok()).collect();
        assert_eq!(locs.len(), 4);
        Some((id, locs[0], locs[1], locs[2], locs[3]))
    };
    let locations: Vec<_> = BufReader::new(File::open(&args[3])?)
        .lines()
        .filter_map(|line| line.ok())
        .filter_map(parse_line)
        .collect();
    let result: Vec<_> = locations
        .par_iter()
        .map(|record| find_nearest_repeat(record, &reference, &alignment))
        .collect();
    for (id, infos) in result {
        let line: Vec<_> = infos.iter().map(|e| format!("{}", e)).collect();
        println!("{}\t{}", id, line.join("\t"));
    }
    Ok(())
}
type Record<'a> = (&'a str, [usize; 7]);
fn find_nearest_repeat<'a>(
    &(ref id, s_a, e_a, s_b, e_b): &'a (String, usize, usize, usize, usize),
    _reference: &[fasta::Record],
    alignments: &[last_tiling::LastTAB],
) -> Record<'a> {
    let nearest_repeat = alignments
        .iter()
        .min_by(|aln1, aln2| {
            let dist1 = dist((s_a, e_a), aln1) + dist((s_b, e_b), aln1);
            let dist2 = dist((s_a, e_a), aln2) + dist((s_b, e_b), aln2);
            dist1.cmp(&dist2)
        })
        .unwrap();
    let info = get_info((s_a, e_a, s_b, e_b), nearest_repeat);
    (id, info)
}

fn dist((s_a, e_a): (usize, usize), aln: &last_tiling::LastTAB) -> usize {
    let ref_len = aln.seq1_len();
    let positions = [
        aln.seq1_start_from_forward(),
        aln.seq1_end_from_forward(),
        aln.seq2_start_from_forward(),
        aln.seq2_end_from_forward(),
    ];
    let distance = |p: usize, q: usize, len: usize| {
        let d = p.max(q) - p.min(q);
        d.min(len - d)
    };
    let min1 = positions
        .iter()
        .map(|p| distance(p, s_a, ref_len))
        .min()
        .unwrap();
    let min2 = positions
        .iter()
        .map(|p| distance(p, e_a, ref_len))
        .min()
        .unwrap();
    min1.min(min2)
}

fn get_info(
    (s_a, e_a, s_b, e_b): (usize, usize, usize, usize),
    aln: &last_tiling::LastTAB,
) -> [usize; 7] {
    [0; 7]
}
