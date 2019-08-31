extern crate bio;
extern crate handmade_bloom_filter;
use handmade_bloom_filter::UpperBoundBFFactory;
use bio::io;
use handmade_bloom_filter::HUGE_MODULO;
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
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    println!("{}",bf);
    Ok(())
}
