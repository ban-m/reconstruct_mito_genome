extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::io::{BufWriter, Write};
use std::path::Path;
const THR: u32 = 500;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut output = BufWriter::new(std::io::stdout());
    for record in bam::Reader::from_path(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .filter(is_doubleclip)
    {
        writeln!(&mut output, "{}", record.cigar().to_string())?;
    }
    Ok(())
}

fn is_doubleclip(record: &bam::Record) -> bool {
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
        start > THR && end > THR
    }
}
