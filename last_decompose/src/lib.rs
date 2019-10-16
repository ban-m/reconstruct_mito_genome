#[macro_use]
extern crate log;
extern crate bio_utils;
extern crate env_logger;
extern crate last_tiling;
extern crate rand;
use bio_utils::fasta;
use last_tiling::LastTAB;
use rand::rngs::StdRng;
use rand::SeedableRng;
mod find_breakpoint;
use find_breakpoint::critical_regions;
use last_tiling::UNIT_SIZE;
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
            None => panic!("{} should be some value, but none.", from),
        };
        let t = match &assignments[to] {
            Some(ref t) => t,
            None => panic!("{} should be some value, but none.", from),
        };
        let new = assignments::merge_two_assignments(f, t);
        assignments[to] = Some(new);
        assignments[from] = None;
    }
    let mut assignments: Vec<_> = assignments.into_iter().filter_map(|e| e).collect();
    let assignment = assignments.pop().unwrap();
    assert!(assignments.is_empty());
    let mut result = vec![vec![]; assignment.get_num_of_cluster()];
    let mut rng: StdRng = SeedableRng::seed_from_u64(2444);
    // The order of read is critical, since the assignments are just array of weight.
    for (idx, read) in read.into_iter().enumerate() {
        let c = assignment.assign(idx, &mut rng);
        result[c].push(read);
    }
    result
}
