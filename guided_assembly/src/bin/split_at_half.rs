extern crate bio;
use bio::io::fasta;
use std::path::Path;
/// This program would output a (multi) fasta file into standard output,
/// where the fasta file consists of each fasta record in the first argument after
/// splitting them at the middle and concatenate the two fragments the other way around.
/// Cleary speking, This function split fasta file like this : <-----A---------C------->
/// into :<-----A----> <-----C------->  and concat like
/// <-----C------->
///            <-----A---->
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let path = &args[1];
    let overlap_len: usize = args[2].parse().unwrap();
    let mut wtr = fasta::Writer::new(std::io::stdout());
    for record in fasta::Reader::from_file(&Path::new(path))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| split_at_half(e, overlap_len))
    {
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok(())
}

fn split_at_half(record: fasta::Record, overlap_len: usize) -> fasta::Record {
    let seq = record.seq();
    let (former, latter) = seq.split_at(seq.len() / 2);
    fasta::Record::with_attrs(
        record.id(),
        record.desc(),
        &[latter, &former[overlap_len..]].concat(),
    )
}
