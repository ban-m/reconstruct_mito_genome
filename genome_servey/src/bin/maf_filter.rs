extern crate bio_utils;
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let threshold: u64 = args[2].parse().unwrap();
    bio_utils::maf::Reader::from_file(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .filter(|record| record.sequence().iter().all(|seq| seq.length() > threshold))
        .for_each(|record| println!("{}", record));
    Ok(())
}
