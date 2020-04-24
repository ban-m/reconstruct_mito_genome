extern crate bio_utils;
extern crate rust_htslib;
use rust_htslib::{bam, bam::Read};
const THR: f64 = 0.09;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mean_cov = mean_coverage(&args[1]);
    let (total, putative) = count_variant_position(&args[1], mean_cov);
    eprintln!("{}\t{}\t{}\t{}", mean_cov, total, putative, THR);
}

fn mean_coverage(file: &str) -> u32 {
    let mut reader = bam::Reader::from_path(file).unwrap();
    let tot = reader
        .pileup()
        .filter_map(|e| e.ok())
        .filter(|p| p.depth() > 3)
        .count();
    let mut reader = bam::Reader::from_path(file).unwrap();
    let sum: u32 = reader
        .pileup()
        .filter_map(|e| e.ok())
        .filter(|p| p.depth() > 3)
        .map(|p| p.depth())
        .sum();
    (sum / tot as u32)
}

fn count_variant_position(file: &str, ave: u32) -> (usize, usize) {
    let mut reader = bam::Reader::from_path(file).unwrap();
    let tot = reader
        .pileup()
        .filter_map(|e| e.ok())
        .filter(|p| p.depth() >= ave / 2)
        .count();
    let mut reader = bam::Reader::from_path(file).unwrap();
    let putative = reader
        .pileup()
        .filter_map(|e| e.ok())
        .filter(|p| p.depth() >= ave / 2)
        .filter(|p| {
            let depth = p.depth();
            let mut count = [0; 4];
            for aln in p.alignments() {
                if !aln.is_del() && !aln.is_refskip() && !aln.record().seq().is_empty() {
                    let base = aln.record().seq()[aln.qpos().unwrap()];
                    let idx = match base {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => continue,
                    };
                    count[idx] += 1;
                }
            }
            let count = minor_allel_count(count);
            count as f64 / depth as f64 > THR
        })
        .count();
    (putative, tot)
}

fn minor_allel_count(mut count: [usize; 4]) -> usize {
    count.sort();
    count[2]
}
