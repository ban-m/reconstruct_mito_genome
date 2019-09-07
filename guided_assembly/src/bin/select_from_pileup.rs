extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::path::Path;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let contig_name = &args[2];
    let start: i32 = args[3].parse().unwrap();
    let end: i32 = args[4].parse().unwrap();
    let mut reader = bam::Reader::from_path(&Path::new(&args[1])).unwrap();
    let header = reader.header();
    let tid = header.tid(contig_name.as_bytes()).unwrap() as i32;
    let mut wtr = bam::Writer::from_stdout(&bam::Header::from_template(header)).unwrap();
    for read in reader
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| is_focal_alignment(e, tid, start, end))
    {
        wtr.write(&read).unwrap();
    }
}

fn is_focal_alignment(align: &bam::Record, tid: i32, start: i32, end: i32) -> bool {
    // Proper alignmnet
    if align.flags() & 0x900 != 0 {
        return false;
    }
    if align.tid() != tid {
        return false;
    }
    let q_start = align.pos() + 100; // safety margin
    let q_end = match align.cigar().end_pos() {
        Ok(res) => res - 100, // same as above
        Err(why) => {
            eprintln!("{:?}", why);
            return false;
        }
    };
    q_start < start && end < q_end
}
