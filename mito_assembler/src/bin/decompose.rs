extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
fn main() -> std::io::Result<()> {
    // read, alignements, contigs, and output path.
    let args: Vec<_> = std::env::args().collect();
    let read = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = bio_utils::fasta::parse_into_vec(&args[3])?;
    let self_alignments = last_tiling::parse_tab_file(&args[4])?;
    let output_path = &args[5];
    let decomposed = last_decompose::decompose(read, alignments, contigs, self_alignments);
    for (idx, reads) in decomposed.into_iter().enumerate() {
        let path = format!("{}/{}.rs", output_path, idx);
        use std::io::BufWriter;
        let mut writer =
            bio_utils::fasta::Writer::new(BufWriter::new(std::fs::File::create(&path)?));
        for r in reads {
            writer.write_record(&r)?;
        }
    }
    Ok(())
}
