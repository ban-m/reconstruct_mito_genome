extern crate bio;
use std::collections::HashSet;
use std::io::{BufReader,BufRead,BufWriter};
use std::path::Path;
use std::fs::File;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let selected:HashSet<String> = BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e|e.ok())
        .filter_map(|line| line.split('\t').nth(0).map(|id|id.to_string()))
        .collect();
    let mut wtr = bio::io::fastq::Writer::new(BufWriter::new(std::io::stdout()));
    for record in bio::io::fastq::Reader::from_file(&Path::new(&args[2]))?
    .records()
        .filter_map(|e|e.ok())
        .filter(|rec| selected.contains(rec.id()))
    {
        wtr.write_record(&record)?;
    }
    Ok(())
}
