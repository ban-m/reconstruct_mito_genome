extern crate env_logger;
extern crate last_decompose;
extern crate log;
extern crate serde;
extern crate serde_json;
use last_decompose::d3_data::convert_to_d3_data_with_assign;
use last_decompose::find_breakpoint::Cluster;
use last_tiling::{Contigs, EncodedRead};
use std::collections::HashMap;
pub mod template;
pub fn dump_viewer(
    results: &HashMap<String, u8>,
    reads: &[EncodedRead],
    clusters: &[Cluster],
    contigs: &Contigs,
) -> std::io::Result<String> {
    let clusters = summarize_clusters(clusters, results);
    let summary = convert_to_d3_data_with_assign(&contigs, &reads, &clusters, &results);
    Ok(serde_json::ser::to_string(&summary).unwrap())
}

fn summarize_clusters(clusters: &[Cluster], results: &HashMap<String, u8>) -> Vec<Cluster> {
    let max = results.values().copied().max().unwrap_or(0) as usize;
    let mut aggregated: Vec<_> = (0..=max)
        .map(|id| Cluster {
            id,
            members: vec![],
        })
        .collect();
    use last_decompose::find_breakpoint::ReadClassify;
    // Put initial clusters into the most probable cluster.
    for cluster in clusters {
        // Collect the reads and clusters
        let reads: std::collections::HashMap<u8, usize> = results
            .iter()
            .filter(|&(ref id, _)| cluster.has(id))
            .map(|e| e.1)
            .fold(HashMap::new(), |mut x, &y| {
                *x.entry(y).or_default() += 1;
                x
            });
        // Determine which clusters these initial cluster should be put on.
        if let Some((&argmax, _)) = reads.iter().max_by_key(|x| x.1) {
            let members = cluster.members.iter().cloned().map(|mut mem| {
                mem.cluster = argmax as usize;
                mem
            });
            aggregated[argmax as usize].members.extend(members);
        }
    }
    aggregated
}
