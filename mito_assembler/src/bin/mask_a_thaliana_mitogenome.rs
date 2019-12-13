const FORWARD_REPEAT: (usize, usize) = (194200, 199400);
const REVERSE_REPEAT: (usize, usize) = (258700, 265300);
extern crate bio_utils;
use bio_utils::fasta;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let record = fasta::parse_into_vec(&args[1])?
        .into_iter()
        .filter(|e| e.id() == "NC_037304.1")
        .nth(0)
        .unwrap();
    let mut seq = record.seq().to_vec();
    let (s, e) = FORWARD_REPEAT;
    for i in s..e {
        seq[i].make_ascii_uppercase();
    }
    let (s, e) = REVERSE_REPEAT;
    for i in s..e {
        seq[i].make_ascii_uppercase();
    }
    let desc = record.desc().map(|e| e.clone());
    let record = fasta::Record::with_data(record.id(), &desc, &seq);
    let mut wtr = fasta::Writer::new(std::io::stdout());
    wtr.write_record(&record)
}
