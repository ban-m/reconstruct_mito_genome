use bio_utils::fasta;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: Vec<_> = fasta::parse_into_vec(&args[1])?;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let clusters: Vec<Vec<usize>> = File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .map(|line| {
            line.split_whitespace()
                .map(|e| e.parse().unwrap())
                .collect()
        })
        .collect();
    for (cluster, ids) in clusters.iter().enumerate() {
        let looped = ids
            .iter()
            .filter(|&&read| {
                reads[read]
                    .desc()
                    .map(|e| e.contains("looped"))
                    .unwrap_or(false)
            })
            .count();
        let reference = ids
            .iter()
            .filter(|&&read| {
                reads[read]
                    .desc()
                    .map(|e| !e.contains("looped"))
                    .unwrap_or(false)
            })
            .count();
        println!("{}\t{}\t{}", cluster, reference, looped);
    }
    Ok(())
}
