extern crate bio;
extern crate handmade_bloom_filter;
extern crate rayon;
use bio::io;
use handmade_bloom_filter::UpperBoundBFFactory;
use handmade_bloom_filter::HUGE_MODULO;
use rayon::prelude::*;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let t: usize = args[2].parse().unwrap();
    let k: usize = args[3].parse().unwrap();
    let input: Vec<_> = io::fastq::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| e.seq().to_vec())
        .collect();
    let bf = UpperBoundBFFactory::default()
        .k(k)
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    eprintln!("{}", bf);
    let input: Vec<_> = io::fastq::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let tiled: Vec<_> = input
        .into_par_iter()
        .map(|record| {
            let (mean, sd) = mean_sd(bf.upper_bound_at_each_position(record.seq()));
            (record.id().to_string(), mean, sd, record.seq().len())
        })
        .collect();
    use std::io::{BufWriter, Write};
    let mut wtr = BufWriter::new(std::io::stdout());
    writeln!(&mut wtr, "ID\tMean\tSD\tLen")?;
    for (id, mean, sd, len) in tiled {
        writeln!(&mut wtr, "{}\t{}\t{}\t{}", id, mean, sd, len)?;
    }
    Ok(())
}

fn mean_sd(xs:Vec<u16>)->(u64,u64){
    let mut sum = 0;
    let mut sumsq = 0;
    let len = xs.len() as u64;
    for x in xs{
        sum += x as u64;
        sumsq += (x as u64).pow(2);
    }
    let ave = sum/len;
    let sd = ((sumsq/len - ave * ave) as f64).sqrt().floor() as u64;
    (ave,sd)
}
