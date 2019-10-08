extern crate last_tiling;
#[macro_use]
extern crate serde;
extern crate serde_json;
const THR: usize = 1_000;
use std::io::BufWriter;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let contig: last_tiling::Contigs = last_tiling::Contigs::from_file(&args[1])?;
    let alns = last_tiling::parse_tab_file(&args[2])?;
    let alns: Vec<_> = alns
        .into_iter()
        .filter(|aln| {
            // Check complete alignment.
            let seq1_cmp = aln.seq1_matchlen() != aln.seq1_len();
            let seq2_cmp = aln.seq2_matchlen() != aln.seq2_len();
            // Check Long alignment.
            let long = aln.seq1_matchlen() > THR && aln.seq2_matchlen() > THR;
            seq1_cmp && seq2_cmp && long
        })
        .map(|aln| {
            let rep1 = Repeat {
                id: contig.get_id(aln.seq1_name()).unwrap(),
                name: aln.seq1_name().to_string(),
                start: aln.seq1_start_from_forward(),
                end: aln.seq1_end_from_forward(),
            };
            let rep2 = Repeat {
                id: contig.get_id(aln.seq2_name()).unwrap(),
                name: aln.seq2_name().to_string(),
                start: aln.seq2_start_from_forward(),
                end: aln.seq2_end_from_forward(),
            };
            RepeatPair { rep1, rep2 }
        })
        .collect();
    let stdout = std::io::stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    serde_json::ser::to_writer(&mut stdout, &alns).unwrap();
    Ok(())
}

#[derive(Serialize, Deserialize)]
struct RepeatPair {
    rep1: Repeat,
    rep2: Repeat,
}

#[derive(Serialize, Deserialize)]
struct Repeat {
    id: u16,
    name: String,
    start: usize,
    end: usize,
}
