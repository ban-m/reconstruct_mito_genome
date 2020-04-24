extern crate bio_utils;
extern crate last_tiling;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let reference = last_tiling::Contigs::from_file(&args[2])?;
    let alignment = last_tiling::parse_tab_file(&args[3])?;
    let encoded_reads = last_tiling::encoding(&reads, &reference, &alignment);
    let unique_units = {
        let mut units: HashMap<_, usize> = HashMap::new();
        for seq in encoded_reads.iter().map(|e| e.seq()) {
            for e in seq.iter().filter_map(|e| e.encode()) {
                *units.entry((e.contig, e.unit)).or_default() += 1;
            }
        }
        units
    };
    let count = unique_units.iter().filter(|x| x.1 > &10).count();
    let mut wtr = bio_utils::fasta::Writer::new(std::fs::File::create(&args[4])?);
    assert_eq!(reads.len(), encoded_reads.len());
    let mut ok_read = 0;
    for (er, read) in encoded_reads.iter().zip(reads) {
        let has_alone_unit = er
            .seq()
            .iter()
            .filter_map(|e| e.encode())
            .any(|e| unique_units[&(e.contig, e.unit)] <= 2);
        if !has_alone_unit {
            ok_read += 1;
            wtr.write_record(&read)?;
        }
    }
    eprintln!("{} => {}", encoded_reads.len(), ok_read);
    println!("{}K", count * last_tiling::UNIT_SIZE / 1_000);
    Ok(())
}
