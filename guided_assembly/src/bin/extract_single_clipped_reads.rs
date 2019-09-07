extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::path::Path;
const THR: u32 = 500;
use std::collections::HashSet;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let record = bam::Reader::from_path(&Path::new(&args[1])).unwrap();
    let header = bam::header::Header::from_template(&record.header());
    let mut output = bam::Writer::from_stdout(&header).unwrap();
    let single_clips: HashSet<Vec<u8>> = bam::Reader::from_path(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .filter(is_singleclip)
        .map(|e| e.qname().to_vec())
        .collect();
    for record in bam::Reader::from_path(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| single_clips.contains(e.qname()))
    {
        output.write(&record).unwrap();
    }
    Ok(())
}

fn is_singleclip(record: &bam::Record) -> bool {
    use rust_htslib::bam::record::Cigar;
    if record.flags() & 0x100 == 0x100 {
        false
    } else if record.flags() & 0x800 == 0x800 {
        false
    } else {
        let (mut start, mut end) = (0, 0);
        for op in record.cigar().iter() {
            match op {
                Cigar::SoftClip(l) | Cigar::HardClip(l) if start == 0 => start += l,
                Cigar::SoftClip(l) | Cigar::HardClip(l) => end += l,
                _ => {}
            }
        }
        start > THR || end > THR
    }
}
