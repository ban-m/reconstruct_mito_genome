extern crate bio;
extern crate bio_utils;
extern crate rand;
extern crate rust_htslib;
extern crate rusty_sandbox;
use rand::{thread_rng};
use rust_htslib::bam::Read;
use std::collections::HashMap;
const THR: u32 = 100;
use rust_htslib::bam::record::Cigar;
pub fn is_contained(record: &rust_htslib::bam::Record) -> bool {
    let cigar = record.cigar();
    let head = match cigar.as_slice().first(){
            Some(Cigar::SoftClip(l)) | Some(Cigar::HardClip(l)) => *l,
        _ => 0,
    };
    let tail = match cigar.as_slice().last(){
            Some(Cigar::SoftClip(l)) | Some(Cigar::HardClip(l)) => *l,
        _ => 0,
    };
    head < THR && tail < THR
}

fn is_overlap(record: &rust_htslib::bam::Record, _ref_length: usize) -> bool {
    let cigar = record.cigar();
    let ref_begin = record.pos() as u32;
    let query_begin = match cigar.as_slice().first() {
        Some(Cigar::SoftClip(l)) | Some(Cigar::HardClip(l)) => *l,
        _ => 0,
    };
    let ref_end = record.cigar().end_pos() as u32;
    let query_end = match cigar.as_slice().last() {
        Some(Cigar::SoftClip(l)) | Some(Cigar::HardClip(l)) => *l,
        _ => 0,
    };
    (ref_begin < THR && query_end < THR) || (query_begin < THR && ref_end < THR)
}
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut sam_reader =
        rust_htslib::bam::Reader::from_path(std::path::Path::new(&args[1])).unwrap();
    let fastq_records: HashMap<_, _> = rusty_sandbox::open_fastq_into_hashmap(&args[2])?;
    let reverse_index = rusty_sandbox::reverse_index(&sam_reader);
    let mut _rng = thread_rng();
    for i in 0..5 {
        let records = sam_reader.records();
        let sam_record = records
            .filter_map(|e| e.ok())
            .filter(|e| e.pos() != 0)
            .filter(|rec| {
                let ref_len = fastq_records[&reverse_index[rec.tid() as usize]]
                    .seq()
                    .len();
                is_overlap(rec, ref_len)
            })
            // .filter(is_contained)
            .nth(i)
            .unwrap();
        let seq1 = &fastq_records[&String::from_utf8_lossy(sam_record.qname()).to_string()];
        let seq2 = &fastq_records[&reverse_index[sam_record.tid() as usize]];
        let seq1 = if sam_record.is_reverse() {
            bio::alphabets::dna::revcomp(seq1.seq())
        } else {
            seq1.seq().to_vec()
        };
        let (seq1_p, ops, seq2_p) = bio_utils::bam::recover_alignment(
            sam_record.cigar().as_slice(),
            &seq1,
            seq2.seq(),
            sam_record.pos() as usize,
        );
        let window = 200;
        for ((line1, line2), line3) in seq1_p
            .chunks(window)
            .zip(ops.chunks(window))
            .zip(seq2_p.chunks(window))
        {
            println!(
                "{}\n{}\n{}\n",
                String::from_utf8_lossy(line1),
                String::from_utf8_lossy(line2),
                String::from_utf8_lossy(line3)
            );
        }
    }
    Ok(())
}
