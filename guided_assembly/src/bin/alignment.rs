extern crate bio;
extern crate rayon;
use bio::alignment::pairwise;
use bio::io::fastq;
use rayon::prelude::*;
use std::path::Path;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let reads: Vec<_> = fastq::Reader::from_file(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let result: Vec<(String, bio::alignment::Alignment, usize)> = reads
        .into_par_iter()
        .map(|read| {
            let mut aligner = pairwise::Aligner::new(-4, -1, |a, b| if a == b { 1 } else { -1 });
            let len = read.seq().len();
            let seq = read.seq();
            let seq2 = bio::alphabets::dna::revcomp(seq);
            let result = aligner.local(seq, &seq2);
            (read.id().to_string(), result, len)
        })
        .collect();
    println!("ID,tstart,tend,rstart,rend,len,score");
    for (id, align, len) in result {
        println!(
            "{},{},{},{},{},{},{}",
            id,
            align.xstart,
            align.xend,
            len - align.yend,
            len - align.ystart,
            align.ylen,
            align.score
        );
    }
    let end = std::time::Instant::now();
    eprintln!("{:?}", end - start);
}
