extern crate bio;
extern crate handmade_bloom_filter;
use handmade_bloom_filter::UpperBoundBFFactory;
use bio::io;
use handmade_bloom_filter::BIG_MODULO;
use std::collections::HashMap;
fn main() -> std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let t:usize = args[2].parse().unwrap();
    let input:Vec<_> = io::fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
    .records()
        .filter_map(|e|e.ok())
        .map(|e|e.seq().to_vec())
        .collect();
    eprintln!("Collect reads:{}",input.len());
    let bf = UpperBoundBFFactory::default().k(20)
        .number_of_hash(5)
        .modulo(BIG_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    let mut hm:HashMap<u16,u32> = HashMap::new();
    for line in &input{
        for count in bf.upper_bound_at_each_position(line){
            hm.entry(count).and_modify(|e| *e +=1).or_insert(1);
        }
    }
    use std::io::{BufWriter,Write};
    let mut wtr = BufWriter::new(std::io::stdout());
    writeln!(&mut wtr,"Count\tNumberOfKmer")?;
    for (count,num_of_kmer) in hm{
        writeln!(&mut wtr, "{}\t{}",count,num_of_kmer)?;
    }
    Ok(())
}
