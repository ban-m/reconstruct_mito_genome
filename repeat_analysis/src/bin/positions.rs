extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
use bio_utils::fasta;
#[macro_use]
extern crate log;
extern crate env_logger;
fn main() -> std::io::Result<()> {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reference = fasta::parse_into_vec(&args[1])?;
    let self_aln = last_tiling::parse_tab_file(&args[2])?;
    let alignments = last_tiling::parse_tab_file(&args[3])?;
    let reads = fasta::parse_into_vec(&args[4])?;
    let contigs = last_tiling::contig::Contigs::new(reference);
    let repeats = last_tiling::into_repeats_full(&self_aln, &contigs);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let encoded_reads: Vec<_> = encoded_reads
        .into_iter()
        .map(last_decompose::ERead::new_no_gapfill)
        .collect();
    let critical_regions = last_decompose::critical_regions(&encoded_reads, &contigs, &repeats);
    use last_decompose::CriticalRegion;
    let num = reads.len() as f64;
    let mut id = 0;
    for c in &critical_regions {
        use last_decompose::find_breakpoint::ReadClassify;
        use last_tiling::UNIT_SIZE;
        if let CriticalRegion::CP(c) = c {
            let count = encoded_reads
                .iter()
                .filter(|read| c.along_with(read))
                .count();
            if count > 10 {
                let (start_a, end_a) = c.contig1().range();
                let (start_a, end_a) = (start_a as usize * UNIT_SIZE, end_a as usize * UNIT_SIZE);
                let (start_b, end_b) = c.contig2().range();
                let (start_b, end_b) = (start_b as usize * UNIT_SIZE, end_b as usize * UNIT_SIZE);
                let frac = count as f64 / num;
                if start_a < start_b {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        id, start_a, end_a, start_b, end_b, count, frac,
                    );
                } else {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        id, start_b, end_b, start_a, end_a, count, frac
                    );
                }
                id += 1;
            }
        }
    }
    Ok(())
}
