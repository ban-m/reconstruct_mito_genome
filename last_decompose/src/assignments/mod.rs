use super::find_breakpoint::*;
use super::last_tiling::{Contigs, EncodedRead};
//use last_tiling::UNIT_SIZE;
use super::find_breakpoint::ReadClassify;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use std::collections::HashSet;
/// if the similatity between two cluster is higher than SIM_THR, an edge would be drawn.
pub const SIM_THR: f64 = 6.0;
/// The number of read spanned.
pub const READ_NUM: usize = 10;
/// The spanning road with reads less than REMOVE_CLUSTER would be discarded.
pub const _REMOVE_CLUSTER: usize = 3;
mod bipartite_matching;
mod dbg_hmms;
mod minimum_spanning_tree;

#[derive(Debug, Clone, Default)]
pub struct Assignment {
    weight: Vec<Vec<f64>>,
    num_of_cluster: usize,
    choises: Vec<usize>,
}

impl Assignment {
    pub fn new(len: usize, num_of_cluster: usize) -> Self {
        let mut weight = vec![vec![0.; len]; num_of_cluster];
        weight[0] = vec![1.; len];
        let choises = (0..num_of_cluster).collect();
        Self {
            weight,
            num_of_cluster,
            choises,
        }
    }
    // TODO
    pub fn push(&mut self, _idx: usize, _weights: &[f64]) {}

    pub fn assign<R: Rng>(&self, idx: usize, r: &mut R) -> usize {
        *self
            .choises
            .choose_weighted(r, |&k| self.weight[k][idx])
            .unwrap()
    }
    pub fn get_num_of_cluster(&self) -> usize {
        self.num_of_cluster
    }
    /// Return the number of READS. Not CLUSTERs.
    /// If you want to get the number of clusters, use `get_num_of_cluster`.
    pub fn len(&self) -> usize {
        self.weight[0].len()
    }
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Calculate similarity between other.
    /// If the readset differ, it would panic.
    fn sim(&self, other: &Self) -> f64 {
        let ka = self.num_of_cluster;
        let kb = self.num_of_cluster;
        self.weight
            .iter()
            .map(|color1| {
                other
                    .weight
                    .iter()
                    .map(|color2| Self::sim_bet_color(color1, color2, ka, kb))
                    .fold(-10_000_000., |max, x| if max < x { x } else { max })
            })
            .fold(-1_000_000., |max, x| if max < x { x } else { max })
    }
    // Return the weight of k-th cluster on each read.
    // In other words, if we let xs be the return value,
    // xs[i] is the weight of cluster k at i-th read.
    fn get_weight_of(&self, k: usize) -> &[f64] {
        &self.weight[k]
    }
    fn sim_bet_color(color1: &[f64], color2: &[f64], ka: usize, kb: usize) -> f64 {
        assert_eq!(color1.len(), color2.len());
        let (ka, kb) = (ka as f64, kb as f64);
        color1
            .iter()
            .zip(color2.iter())
            .map(|(c1, c2)| (ka * c1 - 1.) * (kb * c2 - 1.))
            .fold(0f64, |cum, x| cum + x)
            / ((ka - 1.) * (kb - 1.))
    }
    // Merge self with other, by merging i-th cluster of self
    // to j-th cluster of other, where (i,j) is an element of by.
    // Note that there can be at most one occurence of i and j.
    // In other words, `by` is a matching of self and other.
    // f(i,.)(r) = f_i(r) * k_a / k
    // f(.,j)(r) = f_j(r) * k_b / k
    // f(i,j)(r) = f_i(r) * k_a / k + f_j(r) * k_b / k - 1/k
    // where k_a = # of self's cluster
    // k_b = # of other's cluster, k = # of result cluster.
    fn merge_with_by(&self, other: &Self, by: &[(usize, usize)]) -> Self {
        // Detemine the number of clusters.
        assert_eq!(self.len(), other.len());
        let num_cluster_self = self.get_num_of_cluster();
        let num_cluster_other = other.get_num_of_cluster();
        let num_of_cluster = num_cluster_self + num_cluster_other - by.len();
        let choises: Vec<_> = (0..num_of_cluster).collect();
        let mut has_merged_self = vec![false; num_cluster_self];
        let mut has_merged_other = vec![false; num_cluster_other];
        let mut weight = vec![];
        let k_a = num_cluster_self as f64;
        let k_b = num_cluster_other as f64;
        let k = num_of_cluster as f64;
        for &(i, j) in by {
            has_merged_self[i] = true;
            has_merged_other[j] = true;
            let new_weight: Vec<_> = self
                .get_weight_of(i)
                .iter()
                .zip(other.get_weight_of(j).iter())
                .map(|(&w1, &w2)| (w1 * k_a + w2 * k_b - 1.) / k)
                .collect();
            weight.push(new_weight);
        }
        has_merged_self
            .into_iter()
            .enumerate()
            .filter(|&(_, b)| !b)
            .for_each(|(idx, _)| {
                let w = self
                    .get_weight_of(idx)
                    .iter()
                    .map(|&e| e * k_a / k)
                    .collect::<Vec<_>>();
                weight.push(w);
            });
        has_merged_other
            .into_iter()
            .enumerate()
            .filter(|&(_, b)| !b)
            .for_each(|(idx, _)| {
                let w = other
                    .get_weight_of(idx)
                    .iter()
                    .map(|&e| e * k_b / k)
                    .collect::<Vec<_>>();
                weight.push(w);
            });
        assert_eq!(num_of_cluster, weight.len());
        Assignment {
            weight,
            num_of_cluster,
            choises,
        }
    }
}

// Decomposing the (possibly long) repeat into two. Repeat resolution program.
// Hopefully, it would not be called in the procedure, since if there's no
// mutation inside this repeat, it is no use.
// Also, it is currently able to resolve two repeat. If there are more complex situation,
// It just fails.
// fn local_first_order_decomposing<'a>(
//     reads: &'a [EncodedRead],
//     cr: &ContigPair,
//     contigs: &Contigs,
// ) -> (
//     usize,
//     std::vec::Vec<(usize, usize, &'a EncodedRead)>,
//     std::vec::Vec<(usize, &'a EncodedRead)>,
// ) {
//     let (overlapping_reads, remainings): (Vec<(usize, &EncodedRead)>, Vec<(usize, &EncodedRead)>) =
//         reads
//             .iter()
//             .enumerate()
//             .partition(|(_idx, r)| cr.overlaps_with(r));
//     let (contained_reads, remainings): (Vec<_>, Vec<_>) =
//         remainings.into_iter().partition(|(_, r)| cr.contains(r));
//     // classed_reads[idx] would return the idx-th cluster
//     let (_num_of_cluster, classed_reads): (_, Vec<Vec<_>>) =
//         cr.clustering_reads_into_points(&overlapping_reads);
//     let pairs: Vec<[usize; 2]> = cr.get_competitive_pair(&overlapping_reads);
//     // Let competitive pair be [[i,j],[l,m]]. Thus, we should chose one out of two options:
//     // i connets to l and j connects to m, OR, i connects to m, and j connects to l.
//     // First option is like:
//     // -----i----|->->->|----l----
//     // -----j----|->->->|----m----
//     // Second one is like
//     // -----i----|->->->|----m----
//     // -----j----|->->->|----l----
//     // Note that i-l are 0..4, but we can not tell which to which.
//     let mut rng: StdRng = SeedableRng::seed_from_u64(12123);
//     let separated_reads: Vec<[HashSet<_>; 2]> = pairs
//         .iter()
//         .map(|pair| {
//             let (i, j) = (pair[0], pair[1]);
//             let mut dbg_hmms = dbg_hmms::DeBruijnGraphHiddenMarkovs::new(contigs, 2);
//             classed_reads[i]
//                 .iter()
//                 .for_each(|(_, r)| dbg_hmms.push(0, r));
//             classed_reads[j]
//                 .iter()
//                 .for_each(|(_, r)| dbg_hmms.push(1, r));
//             let (mut class_i, mut class_j) = (HashSet::new(), HashSet::new());
//             contained_reads.iter().for_each(|(idx, read)| {
//                 let weights = dbg_hmms.predict(read);
//                 let assign = *[0, 1].choose_weighted(&mut rng, |&k| weights[k]).unwrap();
//                 if assign == 0 {
//                     class_i.insert(idx);
//                 } else {
//                     class_j.insert(idx);
//                 }
//             });
//             [class_i, class_j]
//         })
//         .collect();
//     // (i->l and j->m) and (i->m and j -> l)
//     let connection_pattern = [[(0, 0, 0), (1, 1, 1)], [(0, 0, 1), (1, 1, 0)]];
//     let selected_pattern = connection_pattern
//         .iter()
//         .max_by_key(|pat| {
//             pat.iter()
//                 .map(|&(x, y, z)| {
//                     separated_reads[x][y]
//                         .intersection(&separated_reads[x][z])
//                         .count()
//                 })
//                 .sum::<usize>()
//         })
//         .unwrap();
//     let classes: Vec<HashSet<usize>> = selected_pattern
//         .iter()
//         .map(|&(x, y, z)| {
//             separated_reads[x][y]
//                 .intersection(&separated_reads[x][z])
//                 .map(|&&e| e)
//                 .collect()
//         })
//         .collect();
//     let mut result: Vec<(usize, usize, _)> = vec![];
//     for (class, y, z) in selected_pattern {
//         result.extend(
//             classed_reads[pairs[0][*y]]
//                 .iter()
//                 .chain(classed_reads[pairs[1][*z]].iter())
//                 .map(|(idx, read)| (*class, *idx, *read)),
//         );
//     }
//     for (idx, read) in contained_reads {
//         if classes[0].contains(&idx) {
//             result.push((0, idx, read));
//         } else if classes[1].contains(&idx) {
//             result.push((1, idx, read));
//         } else {
//             panic!();
//         }
//     }
//     (2, result, remainings)
// }

/// Locally decompose the critical region, and then
/// wave it to remaining reads.
pub fn local_decompose(
    cr: &CriticalRegion,
    reads: &[EncodedRead],
    contigs: &Contigs,
) -> Assignment {
    // Check if there are sufficient number of spannning reads.
    let num_of_spanning_reads = reads.iter().filter(|r| cr.is_spanned_by(r)).count();
    let (num_of_cluster, classed_reads, remaining_reads) = if num_of_spanning_reads > READ_NUM {
        // First, determine the number of clusters
        let (spanning_reads, mut remainings): (
            Vec<_>,
            Vec<_>,
        ) = reads
            .iter()
            .enumerate()
            .partition(|(_idx, r)| cr.is_spanned_by(r));
        // (class, index, read)
        let (num_of_cluster, classed_reads): (usize, Vec<(usize, usize, &EncodedRead)>) =
            cr.separate_reads_into_clusters(spanning_reads);
        remainings.sort_by_key(|(_, r)| cr.distance(r));
        (num_of_cluster, classed_reads, remainings)
    } else {
        // Fallback.
        // Determine each elements.
        // Note: here, the cr should be a instance of ContigPair,
        // as RepeatJunction should have certain amount of spannning reads!
        let cp = match cr {
            CriticalRegion::CP(ref cp) => cp,
            _ => unreachable!(),
        };
        local_first_order_decomposing(reads, cp, contigs)
    };
    // From here, generic classification starts.
    let mut assignment = Assignment::new(reads.len(), num_of_cluster);
    // With default parameter
    let mut dbg_hmms = dbg_hmms::DeBruijnGraphHiddenMarkovs::new(contigs, num_of_cluster);
    // Go on classification.
    for (class, idx, read) in classed_reads {
        let mut weights = vec![0.; num_of_cluster];
        weights[class] = 1.;
        assignment.push(idx, &weights);
        dbg_hmms.push(class, read);
    }
    let mut rng: StdRng = SeedableRng::seed_from_u64(12123);
    for (idx, read) in remaining_reads {
        // Calculate weight
        let weights = dbg_hmms.predict(read);
        assignment.push(idx, &weights);
        // Detemine class.
        let class = *assignment
            .choises
            .choose_weighted(&mut rng, |&k| weights[k])
            .unwrap();
        dbg_hmms.push(class, read);
    }
    assignment
}

/// Return the order in which the assignments should be marged.
/// Note that it could be carried out by first constructing MST,
/// then find center of tree.
/// Maybe we can do more costly procedure?
/// Complexity: O(A^2 * N * K^2), A is number of assignment,
/// N is the number of reads, K is the maximum number of clusters.
/// If we change this proc into iterative minimum merging,
/// the complexity would be O(A^3 * N * K^2).
/// I don't know whether A is large or small, but computation of `sim` would be
/// very costly.
pub fn enumerate_merge_order(assignments: &[Assignment]) -> Vec<(usize, usize)> {
    // First, calculate the similarity between each assignments.
    let mut graph = vec![vec![]; assignments.len()];
    for (i1, a1) in assignments.iter().enumerate() {
        for (i2, a2) in assignments.iter().enumerate() {
            if i1 < i2 {
                let sim = a1.sim(a2);
                // Negativate!
                graph[i1].push((i2, -sim));
                graph[i2].push((i1, -sim));
            }
        }
    }
    // Then, calculate MST.
    let mst = minimum_spanning_tree::mst(graph);
    // Then, find center of MST.
    let center = minimum_spanning_tree::find_center(&mst);
    // Then, dfs and retrieve the edges.
    minimum_spanning_tree::dfs(&mst, center)
}

/// Merge two assignments into a single assignment.
/// Shortly speaking, we first construct a bipatite graph
/// where nodes are clusters in each assignments,
/// and edges with weight w is spanned if the similarity between
/// two clusters are hight enough(SIM_THR>6.0?). SHOULD BE TUNED!
/// Then, we calculate a maximum weight matching, then
/// merge the connected clusters by (f + g) / 2, while
/// let the other assignment as-is.
/// Thus, the final number of clusters would be
/// [number of cluster in as1] + [number of clusters in as2] - [number of clusters merged]
pub fn merge_two_assignments(as1: &Assignment, as2: &Assignment) -> Assignment {
    let nodes_1 = as1.get_num_of_cluster();
    let nodes_2 = as2.get_num_of_cluster();
    // if graph[i] = (j, w),
    // there is a edge BETWEEN node i in as1 and node j in as2 with the weight of w.
    let mut graph = vec![vec![]; nodes_1];
    for i in 0..nodes_1 {
        for j in 0..nodes_2 {
            let w = Assignment::sim_bet_color(
                as1.get_weight_of(i),
                as2.get_weight_of(j),
                nodes_1,
                nodes_2,
            );
            if w > SIM_THR {
                graph[i].push((j, w));
            }
        }
    }
    let edges = bipartite_matching::maximum_weight_matching(nodes_1, nodes_2, &graph);
    as1.merge_with_by(as2, &edges)
}
