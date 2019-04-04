extern crate bio;
extern crate bio_types;
use bio::io::gff;
use bio_types::strand;

fn is_mt(record: &gff::Record) -> bool {
    record.seqname() == "Mt"
}

fn is_gene(record: &gff::Record) -> bool {
    record.feature_type() == "gene"
}

fn to_tab(record: gff::Record) -> (String, u64, u64, char) {
    let start = record.start() - 1;
    let end = *record.end();
    let gene_name = match record.attributes().get("Name") {
        Some(res) => res.to_string(),
        None => record.attributes()["gene_id"].to_string(),
    };
    let strand = match record.strand().unwrap() {
        strand::Strand::Forward => '+',
        strand::Strand::Reverse => '-',
        strand::Strand::Unknown => '-',
    };
    (gene_name, start, end, strand)
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    for (genename, start, end, strand) in
        gff::Reader::from_file(&std::path::Path::new(&args[1]), gff::GffType::GFF3)?
            .records()
            .filter_map(|e| e.ok())
            .filter(is_mt)
            .filter(is_gene)
            .map(to_tab)
    {
        println!(
            "Arabidopsis_thaliana\t{}\t{}\t{}\t{}",
            genename, start, end, strand
        );
    }
    Ok(())
}
