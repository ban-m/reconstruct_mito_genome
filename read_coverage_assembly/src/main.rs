extern crate bio;
use std::io::prelude::*;
use std::path::Path;
use bio::io::fastq;
use std::collections::BTreeSet;
const BUF_CAPACITY:usize= 10_000_000;
const QUANTILE:f64 = 0.75;
fn main() -> std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let ids:BTreeSet<String> = {
        let counts:Vec<(u64,String)> = {
            let mut buf = String::with_capacity(BUF_CAPACITY);
            std::io::stdin().read_to_string(&mut buf)?;
            let mut res:Vec<_> = buf.lines().filter_map(|e|{
                let contents:Vec<&str> = e.split_whitespace().collect();
                let num = contents[0].parse().ok()?;
                Some((num,contents[1].to_string()))
            }).collect();
            res.sort();
            res
        };
        let threshold = counts[(QUANTILE*counts.len() as f64).floor() as usize].0;
        eprintln!("{}",threshold);
        counts.into_iter().filter(|e|e.0 > threshold).map(|e|e.1).collect()
    };
    let mut wtr = fastq::Writer::new(std::io::stdout());
    for record in fastq::Reader::from_file(&Path::new(&args[1]))?
        .records()
        .filter_map(|e|e.ok())
        .filter(|e|ids.contains(e.id())){
            wtr.write_record(&record)?;
        }
    wtr.flush()
}

