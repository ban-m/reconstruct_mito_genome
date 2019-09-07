extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read;
use std::path::Path;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let contig_name = &args[2];
    let start: i32 = args[3].parse().unwrap();
    let mut reader = bam::Reader::from_path(&Path::new(&args[1])).unwrap();
    let header = reader.header();
    let tid = header.tid(contig_name.as_bytes()).unwrap() as i32;
    let length = header.target_len(tid as u32).unwrap() as i32;
    let mut wtr = bam::Writer::from_stdout(&bam::Header::from_template(header)).unwrap();
    for read in reader
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| is_focal_alignment(e, tid, start, length))
    {
        wtr.write(&read).unwrap();
    }
}

fn is_focal_alignment(align: &bam::Record, tid: i32, _start: i32, end: i32) -> bool {
    // Proper alignmnet
    if align.flags() & 0x900 != 0 {
        return false;
    }
    if align.tid() != tid {
        return false;
    }
    if let Some((_head_clip, tail_clip, _q_start, q_end)) = parse_cigar(align) {
        tail_clip < 500 && (end - q_end) < 100 // ending sufficiently near the end.
    } else {
        false
    }
}
// returns Option<(length of head clipping, length of tail clipping,
// the location alignment start, the location alignment ends)>.
fn parse_cigar(align: &bam::Record) -> Option<(u32, u32, i32, i32)> {
    let (mut h, mut m, mut t) = (0, 0, 0);
    for op in align.cigar().iter() {
        match op {
            Cigar::HardClip(l) | Cigar::SoftClip(l) if h == 0 => h += l,
            Cigar::HardClip(l) | Cigar::SoftClip(l) => t += l,
            Cigar::Diff(l)
            | Cigar::Equal(l)
            | Cigar::RefSkip(l)
            | Cigar::Match(l)
            | Cigar::Del(l) => m += l,
            Cigar::Ins(_) | Cigar::Pad(_) => {}
        };
    }
    let start = align.pos();
    Some((h, t, start, start + m as i32))
}
