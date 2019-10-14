extern crate bio_utils;
extern crate last_tiling;
use bio_utils::fasta;
use last_tiling::LastTAB;
use std::collections::HashMap;

use find_breakpoint::critical_regions;

pub fn decompose(
    read: Vec<fasta::Record>,
    alignments: Vec<LastTAB>,
    contigs: Vec<fasta::Record>,
    self_aln:Vec<LastTAB>,
) -> Vec<Vec<fasta::Record>> {
    let contigs = last_tiling::Contigs::new(contigs);
    let repeats = last_tiling::repeat::RepeatPairs::new(&self_alns, &contigs).unwrap();
    let encoded_reads = last_tiling::encoding(&read, &contigs, &alns, &repeats);
    let critical_regions = critical_regions(&encoded_reads, &contigs);
    let mut clusters = vec![];
    for cr in critical_regions {
        clusters.push(local_decompose(cr, &alignments, &read, &contigs));
    }
    let merge_order = enumerate_merge_order();
    for (from, to) in merge_order {
        let new = merge_two_clusters(&clusters[from], &clusters[to]);
        cluster[to] = Some(new);
        cluster[from] = None;
    }
    let clusters = clusters.into_iter().filter_map().collect();
    let cluster = clusters.pop();
    assert!(clusters.is_empty());
    let assignments = vec![vec![]; cluster.len()];
    for (idx, read) in 0..reads.into_iter().enumerate() {
        let c = cluster.assign(idx);
        assignments[c].push(read);
    }
    assignments
}
