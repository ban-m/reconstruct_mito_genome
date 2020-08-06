use bio_utils::fasta;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads: HashMap<String, bool> = fasta::parse_into_vec(&args[1])?
        .into_iter()
        .map(|read| {
            let id = read.id().to_string();
            let is_reference = read.desc().map(|e| !e.contains("looped")).unwrap_or(false);
            (id, is_reference)
        })
        .collect();
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let mut clusters: HashMap<(u64, String), Vec<String>> = HashMap::new();
    for line in File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .skip(1) // Header.
        .filter_map(|e| e.ok())
    {
        let line: Vec<_> = line.split_whitespace().collect();
        let cluster: u64 = line[2].parse().unwrap_or(0);
        let name: String = line[0].to_string();
        let haplotype: String = line[1].to_string();
        clusters.entry((cluster, haplotype)).or_default().push(name);
    }
    for ((cluster_id, haplotype), members) in clusters {
        let looped = members
            .iter()
            .filter_map(|name: &String| reads.get(name))
            .filter(|&&b| !b)
            .count();
        let reference = members
            .iter()
            .filter_map(|name| reads.get(name))
            .filter(|&&b| b)
            .count();
        eprintln!("{}\t{}\t{}\t{}", cluster_id, haplotype, looped, reference);
    }
    Ok(())
}
