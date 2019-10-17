#[macro_use]
extern crate log;
extern crate bio_utils;
extern crate env_logger;
extern crate last_tiling;
extern crate rand;
use log::Level;
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
) -> Vec<Vec<fasta::Record>> {
    let contigs = last_tiling::Contigs::new(contigs);
    // Alignment informations are completely (losslessly) encoded into reads.
    let encoded_reads = last_tiling::encoding(&read, &contigs, &alignments);
    let critical_regions = critical_regions(&encoded_reads, &contigs);
    if log_enabled!(Level::Debug){
        for c in critical_regions{
            debug!("{:?}",c);
        }
        return vec![]
    }
    let assignments: Vec<_> = critical_regions
        .into_iter()
        .map(|cr| assignments::local_decompose(&cr, &encoded_reads, &contigs))
        .collect();
    // We no longer need any annotataion for critical region.
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
