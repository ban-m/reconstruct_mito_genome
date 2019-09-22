extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use env_logger::Env;
fn main() -> std::io::Result<()> {
    // env_logger::from_env(Env::default().default_filter_or("debug")).init();
    env_logger::from_env(Env::default().default_filter_or("warn")).init();
    let args: Vec<_> = std::env::args().collect();
    let peaks: last_tiling::UnitDefinitions =
        serde_json::de::from_reader(std::fs::File::open(&args[1])?)?;
    let reads: Vec<last_tiling::EncodedRead> =
        serde_json::de::from_reader(std::fs::File::open(&args[2])?)?;
    // Check the raw data.
    let reads:Vec<_>  = reads.into_iter().map(|e|e.to_forward().fill_gap())
        .collect();
    Ok(())
}
