extern crate bio;
use bio::io::fasta;
use std::path::Path;
use std::collections::HashMap;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let scaffolds: HashMap<String,_> = fasta::Reader::from_file(&Path::new(&args[1]))
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .map(|e|(e.id().to_string(),e))
        .collect();
    let mut wtr = fasta::Writer::new(std::io::stdout());
    let contig4 = &scaffolds["contig_4"];
    wtr.write(contig4.id(),contig4.desc(),
              &bio::alphabets::dna::revcomp(contig4.seq())).unwrap();
    let contig3 = &scaffolds["contig_3"];
    wtr.write_record(contig3).unwrap();
    for id in vec!["contig_2","contig_5", "contig_6"] {
        wtr.write_record(&scaffolds[id]).unwrap();
    }
}
