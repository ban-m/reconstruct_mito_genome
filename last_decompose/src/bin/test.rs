extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let read = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = bio_utils::fasta::parse_into_vec(&args[3])?;
    let decomposed = last_decompose::decompose(read, alignments, contigs);
    for (idx, reads) in decomposed.into_iter().enumerate() {
        debug!("Cluster {}", idx);
        for r in reads {
            debug!("{}", r.id());
        }
    }
    Ok(())
}
