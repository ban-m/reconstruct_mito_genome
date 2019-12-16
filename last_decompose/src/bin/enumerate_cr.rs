extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use last_decompose::ERead;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = bio_utils::fasta::parse_into_vec(&args[3])?;
    let contigs = last_tiling::Contigs::new(contigs);
    let reads:Vec<_> = last_tiling::encoding(&reads, &contigs, &alignments)
        .into_iter()
        .map(ERead::new)
        .collect();
    let result = last_decompose::critical_regions(&reads, &contigs);
    let wtr = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(wtr.lock());
    serde_json::ser::to_writer_pretty(&mut wtr, &result).unwrap();
    Ok(())
}
