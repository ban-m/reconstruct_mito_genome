extern crate bio;
// use std::path::Path;
// use bio::io;
// use std::fs;
// use std::collections::HashMap;
// use std::io::{BufRead,BufReader};
fn main() -> std::io::Result<()> {
    // let args:Vec<_> = std::env::args().collect();
    // let reads:HashMap<_,_>  = io::fasta::Reader::from_file(&Path::new(&args[1]))?
    // .records()
    //     .filter_map(|e|e.ok())
    //     .filter_map(|e|Some((String::from_utf8(e.id()).ok()?,
    //                          e)))
    //     .collect();
    // let split:Vec<_> = BufReader::new(fs::File::open(&Path::new(&args[2])))
    //     .lines()
    //     .filter_map(|e|e.ok())
    //     .filter_map(|e|{
    //         let e:Vec<_> = e.split('\t').collect();
    //         let pos =
    //     })
    //     ;
    // let mut wtr = io::fasta::Writer::new(std::io::stdout());
    Ok(())
}
