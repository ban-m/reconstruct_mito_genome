extern crate bio;
extern crate rust_htslib;
use bio::io::gff;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::BTreeMap;
use std::io::Result;
use std::path::Path;
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let gff = open_gff(&args[1])?;
    let (mut bam, header) = open_bam(&args[2])?;
    for record in bam
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| e.flags() & 0x900 == 0)
    {
        let tid = record.tid() as u64;
        let refname = &header[&tid];
        if is_in_te(&record, &gff) {
            println!("{}\t{}", refname, 1);
        } else {
            println!("{}\t{}", refname, 0);
        }
    }
    Ok(())
}

fn open_gff(file: &str) -> Result<Vec<gff::Record>> {
    let mut gff: Vec<_> = gff::Reader::from_file(&Path::new(file), gff::GffType::GFF3)?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    gff.sort_by(|a,b| a.start().cmp(&b.start()));
    Ok(gff)
}

fn open_bam(file: &str) -> Result<(bam::Reader, BTreeMap<u64, String>)> {
    let bam = bam::Reader::from_path(&Path::new(file)).unwrap();
    let header = bam.header();
    let names = header.target_names();
    let header: BTreeMap<_, _> = names
        .into_iter()
        .filter_map(|e| {
            Some((
                header.tid(e)? as u64,
                String::from_utf8(e.to_vec()).ok()?,
            ))
        })
        .collect();
    Ok((bam, header))
}

fn is_in_te(record: &bam::Record, gff: &Vec<gff::Record>) -> bool {
    let start = record.pos() as u64;
    let end = record.cigar().end_pos().unwrap() as u64;
    let pos = match gff.binary_search_by(|a| a.start().cmp(&start)) {
        Ok(_) => return true,
        Err(res) => res,
    };
    if pos == 0 {
        gff[pos].start() < &end
    } else if pos == gff.len() {
        start < *gff[pos - 1].end()
    } else {
        start < *gff[pos - 1].end() || gff[pos].start() < &end
    }
}
