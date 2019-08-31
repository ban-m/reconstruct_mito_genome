extern crate bio;
use std::io::{BufReader,BufWriter};
fn main()->std::io::Result<()>{
    // convert input fastq into renamed fastq(encode IDs numerically).
    let mut wtr = bio::io::fastq::Writer::new(BufWriter::new(std::io::stdout()));
    for (idx,record) in bio::io::fastq::Reader::new(BufReader::new(std::io::stdin()))
        .records()
        .filter_map(|e|e.ok())
        .enumerate()
    {
        wtr.write(&format!("{}",idx),record.desc(), record.seq(), record.qual())?;
    }
    Ok(())
}
