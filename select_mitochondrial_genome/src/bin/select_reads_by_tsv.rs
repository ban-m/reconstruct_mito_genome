extern crate bio;
use bio::io::fastq;
use bio::io::fasta;
use std::collections::HashSet;
use std::io::{BufReader,BufRead};
use std::fs::File;
use std::path::Path;


fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let selected_id:HashSet<_> = BufReader::new(File::open(&Path::new(&args[2]))?)
        .lines()
        .filter_map(|e|e.ok())
        .collect();
    let mut stdout = fasta::Writer::new(std::io::stdout());
    for record in fastq::Reader::from_file(&Path::new(&args[1]))?.records().filter_map(|e|e.ok()){
        if selected_id.contains(record.id()){
            stdout.write(record.id(),record.desc(),record.seq())?;
        }
    }
    stdout.flush()?;
    Ok(())
}
