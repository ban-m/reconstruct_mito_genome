extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use env_logger::Env;
fn main() -> std::io::Result<()> {
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    info!("Start");
    let alignments = last_tiling::parse_tab_file(&args[1])?;
    debug!("Alignments:{}", alignments.len());
    let contigs = last_tiling::contig::Contigs::from_file(&args[2])?;
    debug!("Contig files:\n{}", contigs);
    let fasta = bio_utils::fasta::parse_into_vec(&args[3])?;
    let fasta: Vec<_> = fasta
        .into_iter()
        .filter(|e| e.id() == "m54113_160913_184949/12255889/20770_37317")
        .collect();
    let encoded_reads = last_tiling::encoding(&fasta, &contigs, &alignments);
    debug!("Encoded:\t{}", encoded_reads.len());
    for read in encoded_reads.iter() {
        println!("{}", read);
    }
    Ok(())
}
