// extern crate bio;
// extern crate handmade_bloom_filter;
// extern crate rayon;
// use rayon::prelude::*;

// const MODULO: usize = 2_038_081_147; // For a big dataset
// use std::io::{BufWriter, Write};
fn main() -> std::io::Result<()> {
    Ok(())
    // let args: Vec<_> = std::env::args().collect();
    // let path = std::path::Path::new(&args[1]);
    // let k: usize = args[2].parse().unwrap();
    // let lower: usize = args[3].parse().unwrap();
    // if lower < 5 {
    //     panic!("the lower bound should be greater than 5");
    // }
    // let upper: usize = args[4].parse().unwrap();
    // let histo: Vec<usize> = vec![0; upper];
    // for (kmer, count) in count(path, k, lower, upper)? {
    //     histo[count] += 1;
    // }
    // let mut wtr = BufWriter::new(std::io::stdout());
    // for i in lower..upper {
    //     writeln!(&mut wtr, "{}\t{}", i, histo[i])?;
    // }
    // Ok(())
}

// use bio::io::fastq;
// fn count(
//     path: std::path::Path,
//     k: usize,
//     lower: usize,
//     upper: usize,
// ) -> std::io::Result<Vec<(<Vec<u8>, usize>)>> {
//     let input = fastq::Reader::from_file(&path)?;
//     let bf = construct_filter(input,k);
//     let input = fastq::Reader::from_file(&path)?;
//     let records:Vec<_> = get_records(input);
//     Ok(count_inner(input, k, lower, upper))
// }
// use handmade_bloom_filter::BloomFilter;
// fn count_inner<R>(mut reader: fastq::Reader<R>, k:usize, lower:usize, upper:usize)->Vec<(Vec<u8>,usize)>{
//     let mut record = fastq::Record::new();
//     let mut bfs = [BloomFilter::new_with_configures(MODULO, 20, k).unwrap();5];
//     let mut count = HashMap::new();
//     loop{
//         reader.read(&mut record).unwrap();

//     }
// }
