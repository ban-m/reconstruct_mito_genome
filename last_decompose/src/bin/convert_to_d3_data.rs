extern crate env_logger;
extern crate last_decompose;
extern crate last_tiling;
extern crate log;
extern crate serde;
extern crate serde_json;
use std::io::{BufReader, BufWriter};
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let contigs: last_tiling::Contigs =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let reads: Vec<last_tiling::EncodedRead> =
        serde_json::de::from_reader(std::fs::File::open(&args[2]).map(BufReader::new)?).unwrap();
    let repeats: Vec<last_tiling::repeat::RepeatPairs> = last_tiling::repeat::open(&args[3])?;
    let alns = last_tiling::parse_tab_file(&args[4])?;
    let clusters = {
        let reads: Vec<_> = reads
            .iter()
            .map(|read| last_decompose::ERead::new(read.clone()))
            .collect();
        last_decompose::initial_clusters(&reads, &contigs, &repeats, &alns)
    };
    let summary = last_decompose::d3_data::convert_to_d3_data(&contigs, &reads, &clusters);
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    serde_json::ser::to_writer(&mut stdout, &summary).unwrap();
    Ok(())
}
