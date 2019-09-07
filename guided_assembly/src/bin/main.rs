extern crate bio;
use bio::io::fastq;
use std::collections::HashSet;
use std::io::Read;
use std::path::Path;
use std::time;
const BUFFER_SIZE: usize = 10_000_000_000;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let now = time::Instant::now();
    let paf: HashSet<String> = {
        let mut buffer = String::with_capacity(BUFFER_SIZE);
        let mut input = std::io::stdin();
        input.read_to_string(&mut buffer)?;
        buffer.lines().map(|e| e.to_string()).collect()
    };
    let dur = time::Instant::now() - now;
    eprintln!("Read Names :{},{}", dur.as_secs(), dur.subsec_nanos());
    let mut wtr = fastq::Writer::new(std::io::stdout());
    for read in fastq::Reader::from_file(&Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
    {
        if paf.contains(read.id()) {
            wtr.write_record(&read)?;
        }
    }
    wtr.flush()
}
