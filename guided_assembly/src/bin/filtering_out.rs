extern crate rust_htslib;
use std::collections::HashSet;
use rust_htslib::bam;
use rust_htslib::bam::Read;
const FRACTION: i32 = 10;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut input = bam::Reader::from_path(&std::path::Path::new(&args[1])).unwrap();
    let header = input.header();
    let mut wtr = bam::Writer::from_path(
        &std::path::Path::new(&args[2]),
        &bam::Header::from_template(&header),
    )
    .unwrap();
    for record in input
        .records()
        .filter_map(|e| e.ok())
        .filter(good_alignment)
    {
        wtr.write(&record).unwrap();
    }
}

fn good_alignment(record: &bam::Record) -> bool {
    if record.mapq() < 50 || record.flags() & 0x900 != 0 {
        return false;
    }
    let (start, end, length) = get_start_end_length(record);
    if (end - start).abs() < length / FRACTION {
        return false;
    }
    true
}
use rust_htslib::bam::record::Cigar;
fn get_start_end_length(record: &bam::Record) -> (i32, i32, i32) {
    let (start, mut end, mut length) = (record.pos(), record.pos(), 0);
    for op in record.cigar().iter() {
        match op {
            Cigar::Ins(l) | Cigar::Diff(l) | Cigar::Match(l) | Cigar::Equal(l) => {
                end += *l as i32;
                length += *l as i32;
            }
            Cigar::HardClip(l) | Cigar::SoftClip(l) => length += *l as i32,
            Cigar::Del(_) | Cigar::Pad(_) | Cigar::RefSkip(_) => (),
        }
    }
    (start, end, length)
}
