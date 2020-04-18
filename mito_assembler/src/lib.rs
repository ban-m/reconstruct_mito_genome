extern crate env_logger;
extern crate last_decompose;
#[macro_use]
extern crate log;
extern crate serde;
extern crate serde_json;
use last_decompose::d3_data::convert_to_d3_data_with_assign;
use last_decompose::find_breakpoint::Cluster;
use last_decompose::find_breakpoint::COVERAGE_THR;
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
    for cl in clusters.iter() {
        debug!("{:?}", cl);
    }
    let summary = convert_to_d3_data_with_assign(&contigs, &reads, &clusters, &results);
    Ok(serde_json::ser::to_string(&summary).unwrap())
}

fn summarize_clusters(clusters: &[Cluster], results: &HashMap<String, u8>) -> Vec<Cluster> {
    let max = results.values().copied().max().unwrap_or(0) as usize;
    (0..=max)
        .map(|cl| {
            let reads: std::collections::HashSet<_> = results
                .iter()
                .filter(|&(_, &e)| e == cl as u8)
                .map(|e| e.0)
                .cloned()
                .collect();
            let members: Vec<_> = clusters
                .iter()
                .filter(|cl| cl.ids().intersection(&reads).count() > COVERAGE_THR)
                .flat_map(|cluster| {
                    cluster.members.iter().map(|e| {
                        use last_decompose::find_breakpoint::Member;
                        Member {
                            cr: e.cr.clone(),
                            cluster: cl,
                        }
                    })
                })
                .collect();
            let id = cl;
            Cluster { id, members }
        })
        .collect()
}
