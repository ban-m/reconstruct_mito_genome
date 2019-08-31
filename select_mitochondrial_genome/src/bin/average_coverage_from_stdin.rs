extern crate rayon;
extern crate select_mitochondrial_genome;
use rayon::prelude::*;
use select_mitochondrial_genome::ReadSummary;
use std::io::Result;

fn main() -> Result<()> {
    let summary_of_each_read = summarize_coverage();
    use std::io::{BufWriter, Write};
    let mut wtr = BufWriter::new(std::io::stdout());
    for (id, summary) in summary_of_each_read {
        writeln!(
            &mut wtr,
            "{}\t{}\t{}\t{}",
            id, summary.mean, summary.sd, summary.length
        )?;
    }
    Ok(())
}

fn summarize_coverage() -> Vec<(String, ReadSummary)> {
    select_mitochondrial_genome::stdin_to_intervals()
        .into_par_iter()
        .map(|(id, interval)| (id, ReadSummary::from_interval(interval)))
        .collect()
}
