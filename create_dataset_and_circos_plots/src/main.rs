extern crate bio_utils;
use bio_utils::fasta;
const SPLIT_AT: usize = 10_000;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let record = &fasta::parse_into_vec(&args[1]).unwrap()[0];
    let (head, tail) = record.seq().split_at(SPLIT_AT);
    let seq: Vec<_> = tail.iter().chain(head.iter()).copied().collect();
    let record = fasta::Record::with_data(
        record.id(),
        &Some("Split at 10_000 position".to_string()),
        &seq,
    );
    let mut wtr = fasta::Writer::new(std::io::stdout());
    wtr.write_record(&record).unwrap();
}
