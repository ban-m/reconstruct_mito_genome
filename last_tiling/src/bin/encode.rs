extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use env_logger::Env;
use std::io::Write;
fn main() -> std::io::Result<()> {
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    info!("Start");
    let alignments = last_tiling::parse_tab_file(&args[1])?;
    info!("Alignments:{}", alignments.len());
    let contigs = last_tiling::contig::Contigs::from_file(&args[2])?;
    info!("Contig files:\n{}", contigs);
    let fasta = bio_utils::fasta::parse_into_vec(&args[3])?;
    info!("Read num\t{}", fasta.len());
    let repeats = last_tiling::repeat::open(&args[4])?;
    info!("Repeats:{:?}", repeats.len());
    // let alignments = last_tiling::remove_repeats(alignments, &contigs, &repeats);
    // info!("Filter repeat:{}", alignments.len());
    let encoded_reads = last_tiling::encoding_w_repeat(&fasta, &contigs, &alignments, &repeats);
    info!("Encoded:\t{}", encoded_reads.len());
    let mut wtr = std::fs::File::create(&args[5])?;
    wtr.write_all(serde_json::ser::to_string_pretty(&contigs)?.as_bytes())?;
    let mut wtr = std::fs::File::create(&args[6])?;
    wtr.write_all(serde_json::ser::to_string(&encoded_reads)?.as_bytes())?;
    Ok(())
}
