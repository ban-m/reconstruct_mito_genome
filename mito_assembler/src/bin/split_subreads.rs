// A program to split subreads from a fasta file from PacBio sequencer.
// This is needed to make the sample rate be uniform.
extern crate bio_utils;
use bio_utils::fasta;
fn main() -> std::io::Result<()> {
    let stdin = std::io::stdin();
    pacbio_fasta(fasta::parse_into_vec_from(stdin.lock())?)
}

fn pacbio_fasta(records: Vec<fasta::Record>) -> std::io::Result<()> {
    let stdout = std::io::stdout();
    let mut stdout = fasta::Writer::new(stdout.lock());
    let (mut tot, mut rem) = (0, 0);
    let mut current_well = String::new();
    let mut buffer: Vec<_> = Vec::with_capacity(10);
    let parse_well = |id: &str| match id.rsplitn(2, '/').nth(1) {
        Some(res) => res.to_string(),
        None => {
            eprintln!("{}", id);
            "".to_string()
        }
    };
    let mut flush = |buffer: &mut Vec<fasta::Record>| -> () {
        if let Some(res) = buffer.into_iter().max_by_key(|e| e.seq().len()) {
            rem += res.seq().len();
            stdout.write_record(&res).unwrap();
        }
        buffer.clear();
    };
    for read in records {
        tot += read.seq().len();
        let well = parse_well(read.id());
        if well != current_well {
            // Update current well
            current_well = well;
            // Output
            flush(&mut buffer);
        }
        buffer.push(read);
    }
    flush(&mut buffer);
    eprintln!(
        "Collected. All:{}Gb Remaining:{}Gb",
        tot / 1_000_000_000,
        rem / 1_000_000_000
    );
    Ok(())
}
