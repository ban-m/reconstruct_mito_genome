use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy)]
struct Edge {
    from: usize,
    to: usize,
    weight: f64,
}

use std::cmp;
use std::cmp::Ordering;
impl cmp::Ord for Edge {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.weight.partial_cmp(&other.weight) {
            Some(Ordering::Equal) => Ordering::Equal,
            Some(Ordering::Less) => Ordering::Greater,
            Some(Ordering::Greater) => Ordering::Less,
            _ => panic!(),
        }
    }
}

impl cmp::PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl cmp::Eq for Edge {}

impl cmp::PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
    }
}

impl Edge {
    fn new(from: usize, to: usize, weight: f64) -> Self {
        Self { from, to, weight }
    }
}

pub fn mst(graph: Vec<Vec<(usize, f64)>>) -> Vec<Vec<(usize, f64)>> {
    let size = graph.len();
    let mut is_in_tree: Vec<bool> = vec![false; size];
    let mut num_of_nodes = 0;
    let mut edges = BinaryHeap::new();
    let mut min_weight = 0.;
    let mut used_edges = vec![vec![]; size];
    // Initialize
    is_in_tree[0] = true;
    num_of_nodes += 1;
    for &(to, weight) in &graph[0] {
        edges.push(Edge::new(0, to, weight));
    }
    while num_of_nodes < size {
        assert!(min_weight <= 0.);
        let Edge { from, to, weight } = edges.pop().unwrap();
        if is_in_tree[to] {
            continue;
        } else {
            // use this edge.
            is_in_tree[to] = true;
            num_of_nodes += 1;
            min_weight += weight;
            for &(to, weight) in &graph[to] {
                edges.push(Edge::new(from, to, weight));
            }
            used_edges[from].push((to, weight));
            used_edges[to].push((from, weight));
        }
    }
    used_edges
}

// Very small constant.
const SMALL: f64 = -100_000_000_000.0;

/// Find the center of MST. In other words,
/// argmin_i max_j D[i,j], where D[i,j] is
/// the shortest distance from i to j.
/// Note that here, the distance is always negative.
/// It is very wierd, but works.
pub fn find_center(mst: &[Vec<(usize, f64)>]) -> usize {
    let (x, _) =
        far_reaching(&mst, 0)
            .into_iter()
            .enumerate()
            .fold(
                (0, SMALL),
                |(idx, max), (i, x)| if max < x { (i, x) } else { (idx, max) },
            );
    let x_dfs = far_reaching(&mst, x);
    let (y, _) = x_dfs
        .clone()
        .into_iter()
        .enumerate()
        .fold(
            (0, SMALL),
            |(idx, max), (i, x)| if max < x { (i, x) } else { (idx, max) },
        );
    let y_dfs = far_reaching(&mst, y);
    // The longest distance from i is
    // either the distance from x or from y,
    // where x and y are the most distant pair.
    // The distance should be negative.
    (0..mst.len())
        .enumerate()
        .map(|(idx, i)| {
            let (x, y) = (x_dfs[i], y_dfs[i]);
            if y < x {
                (idx, x)
            } else {
                (idx, y)
            }
        })
        .fold(
            (0, 1_000_000_000.),
            |(idx, min), (i, x)| if x < min { (i, x) } else { (idx, min) },
        )
        .0
}

// Calculate vector xs where
// xs[i] = the maximum distance from start to i-th node.
fn far_reaching(graph: &[Vec<(usize, f64)>], start: usize) -> Vec<f64> {
    let size = graph.len();
    let mut is_reached: Vec<bool> = vec![false; size];
    // Anything is ok.
    let mut dist = vec![-1.; size];
    let mut stack = vec![start];
    dist[start] = 0.0;
    'dfs: while !stack.is_empty() {
        let from = *stack.last().unwrap();
        if !is_reached[from] {
            is_reached[from] = true;
        }
        for &(to, weight) in &graph[from] {
            if !is_reached[to] {
                // It can be updated with no comparison(it's tree!).
                dist[to] = dist[from] + weight;
                stack.push(to);
                continue 'dfs;
            }
            // For each edgse.
        }
        // Postorder
        let _from = stack.pop().unwrap();
    }
    dist.iter().for_each(|&e| assert!(e <= 0.));
    dist
}

/// Depth first search and retuen the order by which nodes can be shurinked.
/// Each element of return values should be (from,to), where
/// we can merge `from`-th node into `to`-th node.
/// Thus, the edges and orders are ALREADY reversed.
pub fn dfs(graph: &[Vec<(usize, f64)>], start: usize) -> Vec<(usize, usize)> {
    let size = graph.len();
    let mut is_reached = vec![false; size];
    let mut order = vec![];
    let mut stack = vec![start];
    'dfs: while !stack.is_empty() {
        let parent = *stack.last().unwrap();
        if !is_reached[parent] {
            is_reached[parent] = true;
        }
        for &(child, _) in &graph[parent] {
            if !is_reached[child] {
                // Margable from child to parent.
                order.push((child, parent));
                stack.push(child);
                continue 'dfs;
            }
            // For each edge
        }
        // Postorder
        let _parent = stack.pop().unwrap();
    }
    order.reverse();
    order
}
