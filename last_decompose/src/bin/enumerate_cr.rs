extern crate bio_utils;
extern crate env_logger;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate serde;
extern crate serde_json;
use last_decompose::ERead;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = last_tiling::Contigs::from_file(&args[3])?;
    let reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let reads: Vec<_> = reads.into_iter().map(ERead::new_no_gapfill).collect();
    let res = last_decompose::find_breakpoint::confluent_position(&alignments, &contigs, 150);
    for e in res {
        debug!("{:?}", e);
    }
    let res = last_decompose::find_breakpoint::contigpair_position(&reads, &contigs);
    for e in res {
        debug!("{:?}", e);
    }
    Ok(())
}
