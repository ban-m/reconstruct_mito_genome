extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::path::Path;
use std::cmp::Ordering;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut bam = bam::Reader::from_path(&Path::new(&args[1])).unwrap();
    let header = bam.header();
    let names = hash_names(header);
    let mut result:Vec<_> = bam.pileup()
        .filter_map(|e| e.ok())
        .map(|e| to_start_and_stop(&e,&names))
        .collect();
    result.sort_by(|a,b|
        match (a.0).cmp(b.0){
            Ordering::Equal => (a.1).cmp(&b.1),
            x => x,
        });
    println!("refname\tposition\tdepth\tnumber_of_start_read\tnumber_of_stop_read");
    for (refname, pos, depth, start, stop) in result{
        println!(
            "{}\t{}\t{}\t{}\t{}",
            refname, pos, depth, start, stop
            );
    }
    Ok(())
}

fn hash_names(header:&bam::HeaderView)->HashMap<u32,String>{
    let names = header.target_names();
    names.into_iter()
        .filter_map(|name|{
            let tid = header.tid(name)?;
            let name = String::from_utf8_lossy(name).to_string();
            Some((tid, name))
        })
        .collect()
}
fn to_start_and_stop<'a>(pileup: &bam::pileup::Pileup, header:&'a HashMap<u32,String>) -> (&'a str, usize, u32, u32, u32){
    let refname = &header[&pileup.tid()];
    let pos = pileup.pos() as usize;
    let to_int = |x| if x {1}else{0};
    let depth = pileup.depth();
    let (start,end) = pileup.alignments()
        .map(|e|(to_int(e.is_head()),to_int(e.is_tail())))
        .fold((0,0),|(s,e),(t,f)|(s+t,e+f));
    (refname, pos, depth, start,end)
}
