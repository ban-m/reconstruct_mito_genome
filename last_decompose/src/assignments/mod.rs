use super::find_breakpoint::BroadRepeat;
use super::last_tiling::{repeat::RepeatPairs, Contigs, EncodedRead, LastTAB};
use rand::seq::SliceRandom;
use rand::Rng;
/// if the similatity between two cluster is higher than SIM_THR, an edge would be drawn.
pub const SIM_THR: f64 = 6.0;
mod bipartite_matching;
mod minimum_spanning_tree;
#[derive(Debug, Clone, Default)]
pub struct Assignment {
    weight: Vec<Vec<f64>>,
    num_of_cluster: usize,
    choises: Vec<usize>,
}

impl Assignment {
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
                    .fold(-10000000., |max, x| if max < x { x } else { max })
            })
            .fold(-1000000., |max, x| if max < x { x } else { max })
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
        for &(i, j) in by {
            has_merged_self[i] = true;
            has_merged_other[j] = true;
            let new_weight: Vec<_> = self
                .get_weight_of(i)
                .iter()
                .zip(other.get_weight_of(j).iter())
                .map(|(&w1, &w2)| (w1 + w2) / 2.)
                .collect();
            weight.push(new_weight);
        }
        has_merged_self
            .into_iter()
            .enumerate()
            .filter(|&(_, b)| !b)
            .for_each(|(idx, _)| weight.push(self.get_weight_of(idx).to_vec()));
        has_merged_other
            .into_iter()
            .enumerate()
            .filter(|&(_, b)| !b)
            .for_each(|(idx, _)| weight.push(other.get_weight_of(idx).to_vec()));
        assert_eq!(num_of_cluster, weight.len());
        Assignment {
            weight,
            num_of_cluster,
            choises,
        }
    }
}

#[allow(unused_variables)]
pub fn local_decompose(
    cr: BroadRepeat,
    alns: &[LastTAB],
    reads: &[EncodedRead],
    contigs: &Contigs,
    repeat: &[RepeatPairs],
) -> Assignment {
    Assignment::default()
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
