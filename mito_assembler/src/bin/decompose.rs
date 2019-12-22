extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
fn main() -> std::io::Result<()> {
    // read, alignements, contigs, and output path.
    let args: Vec<_> = std::env::args().collect();
    let read = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = bio_utils::fasta::parse_into_vec(&args[3])?;
    // Deprecated options.
    let _self_alignments = last_tiling::parse_tab_file(&args[4])?;
    let output_path = &args[5];
    let decomposed = last_decompose::decompose(read, alignments, contigs, vec![]);
    for (_idx, _reads) in decomposed.into_iter().enumerate() {
        let _path = format!("{}/{}.rs", output_path, _idx);
    }
    Ok(())
}
