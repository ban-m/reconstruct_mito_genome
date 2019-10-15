extern crate bio_utils;
extern crate last_tiling;
extern crate rand;
use bio_utils::fasta;
use last_tiling::LastTAB;
use rand::SeedableRng;
use rand::rngs::StdRng;
mod find_breakpoint;
use find_breakpoint::critical_regions;

mod assignments;

pub fn decompose(
    read: Vec<fasta::Record>,
    alignments: Vec<LastTAB>,
    contigs: Vec<fasta::Record>,
    self_alns: Vec<LastTAB>,
) -> Vec<Vec<fasta::Record>> {
    let contigs = last_tiling::Contigs::new(contigs);
    let repeats = last_tiling::into_repeats(&self_alns, &contigs);
    let encoded_reads = last_tiling::encoding(&read, &contigs, &alignments, &repeats);
    let critical_regions = critical_regions(&encoded_reads, &contigs);
    let mut assignments = vec![];
    for cr in critical_regions {
        assignments.push(assignments::local_decompose(
            cr,
            &alignments,
            &encoded_reads,
            &contigs,
            &repeats,
        ));
    }
    let merge_order = assignments::enumerate_merge_order(&assignments);
    let mut assignments: Vec<_> = assignments.into_iter().map(|e| Some(e)).collect();
    for (from, to) in merge_order {
        // Merge `from` assignment into `to`.
        let f = match &assignments[from] {
            Some(ref f) => f,
            None => continue,
        };
        let t = match &assignments[to] {
            Some(ref t) => t,
            None => continue,
        };
        let new = assignments::merge_two_assignments(f, t);
        assignments[to] = Some(new);
        assignments[from] = None;
    }
    let mut assignments:Vec<_> = assignments.into_iter().filter_map(|e| e).collect();
    let assignment = assignments.pop().unwrap();
    assert!(assignments.is_empty());
    let mut result = vec![vec![]; assignment.len()];
    let mut rng: StdRng = SeedableRng::seed_from_u64(2444);
    for (idx, read) in read.into_iter().enumerate() {
        let c = assignment.assign(idx,&mut rng);
        result[c].push(read);
    }
    result
}
