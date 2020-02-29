extern crate bio_utils;
extern crate last_tiling;
extern crate rayon;
use bio_utils::fasta;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reference = fasta::parse_into_vec(&args[1])?;
    let alignment: Vec<_> = last_tiling::parse_tab_file(&args[2])?;
    let parse_line = |line: String| {
        let line: Vec<&str> = line.split('\t').collect();
        let id = line[0].to_string();
        let locs: Vec<usize> = line[1..6].iter().filter_map(|e| e.parse().ok()).collect();
        assert_eq!(locs.len(), 5);
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
            let dist1 = dist((s_a, e_a, s_b, e_b), aln1);
            let dist2 = dist((s_a, e_a, s_b, e_b), aln2);
            dist1.cmp(&dist2)
        })
        .unwrap();
    let info = get_info((s_a, e_a, s_b, e_b), nearest_repeat);
    (id, info)
}

fn distance((x, y): (usize, usize), (p, q): (usize, usize), len: usize) -> usize {
    // x < y and p < q
    // Overlapping detection
    if y < p {
        // [x,y] ... [p,q]
        (p - y).min(len - (q - x))
    } else if q < x {
        // [p, q] ... [x,y]
        (x - q).min(len - (y - p))
    } else {
        // Overlapping
        0
    }
}

fn dist((s_a, e_a, s_b, e_b): (usize, usize, usize, usize), aln: &last_tiling::LastTAB) -> usize {
    let ref_len = aln.seq1_len();
    let pos = [
        (aln.seq1_start_from_forward(), aln.seq1_end_from_forward()),
        (aln.seq2_start_from_forward(), aln.seq2_end_from_forward()),
    ];
    let min1 = distance(pos[0], (s_a, e_a), ref_len) + distance(pos[1], (s_b, e_b), ref_len);
    let min2 = distance(pos[1], (s_a, e_a), ref_len) + distance(pos[0], (s_b, e_b), ref_len);
    min1.min(min2)
}

fn get_info(
    (s_a, e_a, s_b, e_b): (usize, usize, usize, usize),
    aln: &last_tiling::LastTAB,
) -> [usize; 7] {
    let mut infos = [0; 7];
    // Determine pairing.
    let ref_len = aln.seq1_len();
    let pos = [
        (aln.seq1_start_from_forward(), aln.seq1_end_from_forward()),
        (aln.seq2_start_from_forward(), aln.seq2_end_from_forward()),
    ];
    let min1 = distance(pos[0], (s_a, e_a), ref_len) + distance(pos[1], (s_b, e_b), ref_len);
    let min2 = distance(pos[1], (s_a, e_a), ref_len) + distance(pos[0], (s_b, e_b), ref_len);
    if min1 < min2 {
        // A-> alignment1, B-> alignment2
        infos[0] = aln.seq1_start_from_forward();
        infos[1] = aln.seq1_end_from_forward();
        infos[2] = aln.seq2_start_from_forward();
        infos[3] = aln.seq2_end_from_forward();
    } else {
        // A-> alignment2, B-> alignment1
        infos[0] = aln.seq2_start_from_forward();
        infos[1] = aln.seq2_end_from_forward();
        infos[2] = aln.seq1_start_from_forward();
        infos[3] = aln.seq1_end_from_forward();
    }
    infos[4] = aln.seq1_matchlen() + aln.seq2_matchlen();
    infos[5] = min1;
    infos[6] = min2;
    infos
}
