extern crate last_decompose;
extern crate serde_json;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let alignments: Vec<_> = last_tiling::parse_tab_file(&args[1])?;
    let reference = bio_utils::fasta::parse_into_vec(&args[2])?;
    use last_decompose::annotate_contigs_to_reference::annotate_aln_contigs_to_ref;
    let summary = annotate_aln_contigs_to_ref(&reference, &alignments);
    let result = serde_json::ser::to_string(&summary).unwrap();
    println!("{}", result);
    Ok(())
}
