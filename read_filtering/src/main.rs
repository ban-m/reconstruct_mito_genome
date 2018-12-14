extern crate bio;
use bio::io;
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args[1] == "pacbio" {
        if args[2].ends_with('a') {
            pacbio_fasta(&args[2])?;
        } else if args[2].ends_with('q') {
            pacbio_fastq(&args[2])?;
        } else {
            eprintln!("Please input fasta or fastq file as a second argument. It should be end 'a' or 'q'")
        }
    } else if args[1] == "nanopore" {
        eprintln!("Sorry, currently the procedure for nanopore has not been implemented");
    } else {
        eprintln!("Please give the sequencer name as a first argument: nanopore or pacbio.");
    }
    Ok(())
}

fn pacbio_fasta(fasta: &str) -> std::io::Result<()> {
    let (mut tot,mut rem) = (0,0);
    let mut current_well = String::new();
    let mut buffer :Vec<_> = Vec::with_capacity(10);
    let mut stdout = io::fasta::Writer::new(std::io::stdout());
    let parse_well = |id:&str|{
        match id.rsplitn(2,'-').nth(1){
            Some(res) => res.to_string(),
            None => {eprintln!("{}",id);"".to_string()},
        }
    };
    let mut flush = |buffer:&mut Vec<io::fasta::Record>|->() {
        if let Some(res) = buffer.into_iter().max_by_key(|e|e.seq().len()){
            rem += res.seq().len();
            stdout.write_record(&res).unwrap();
        }
        buffer.clear();
    };
    for read in io::fasta::Reader::from_file(&Path::new(fasta))?
    .records()
        .filter_map(|e| e.ok()){
        tot += read.seq().len();
            let well = parse_well(read.id());
            if well != current_well{
                // Update current well
                current_well = well;
                // Output
                flush(&mut buffer);
            }
            buffer.push(read);
        }
    flush(&mut buffer);
    stdout.flush().unwrap();
    eprintln!("Collected. All:{}Gb Remaining:{}Gb",tot/1_000_000_000,rem/1_000_000_000);
    Ok(())
}

fn pacbio_fastq(fastq: &str) -> std::io::Result<()> {
    let (mut tot,mut rem) = (0,0);
    let mut current_well = String::new();
    let mut buffer :Vec<_> = Vec::with_capacity(10);
    let mut stdout = io::fastq::Writer::new(std::io::stdout());
    let parse_well = |id:&str|{
        match id.rsplitn(2,'/').nth(1){
            Some(res) => res.to_string(),
            None => {eprintln!("{}",id);"".to_string()},
        }
    };
    let mut flush = |buffer:&mut Vec<io::fastq::Record>|->() {
        if let Some(res) = buffer.into_iter().max_by_key(|e|e.seq().len()){
            rem += res.seq().len();
            stdout.write_record(&res).unwrap();
        }
        buffer.clear();
    };
    for read in  io::fastq::Reader::from_file(&Path::new(fastq))?
    .records()
        .filter_map(|e| e.ok()){
            tot += read.seq().len();
            let well = parse_well(read.id());
            eprintln!("{}",well);
            if well != current_well{
                // Update current well
                current_well = well;
                // Output
                flush(&mut buffer);
            }
            buffer.push(read);
        }
    flush(&mut buffer);
    stdout.flush().unwrap();
    eprintln!("Collected. All:{}Gb Remaining:{}Gb",tot/1_000_000_000,rem/1_000_000_000);
    Ok(())
}
