// Find maximum weight matching between the given bipartite matching.
// First, it converts the given graph into flow graph,
// where the const is the negative weight of each edges,
// and capacity is just 1.
// Then, the minimum cost flow is corresponding to the maximum matching of the original graph.
// Also, the flow at each edges either 1 or 0(can be proven).
// Note that the `original_graph` should be a adjacency list. In other words,
// If there is only one edge from index 1 of nodes_1 to index 2 of nodes_2,
// then original_graph[1] should be vec![(2,10.)] or something like that.
pub fn maximum_weight_matching(
    nodes_1: usize,
    nodes_2: usize,
    original_graph: &[Vec<(usize, f64)>],
) -> Vec<(usize, usize)> {
    let total_nodes = nodes_1 + nodes_2 + 2; // for start and end node.
    let mut graph = vec![vec![]; total_nodes];
    let start = 0;
    let end = nodes_1 + nodes_2 + 1;
    for (i, edges) in original_graph.iter().enumerate() {
        for (j, weight) in edges {
            // Make the edges negative.
            graph[i + 1].push((j + nodes_1 + 1, 1, -weight));
        }
    }
    // Add edges from 0(start)-> nodes in nodes_1
    for i in 1..nodes_1 + 1 {
        graph[0].push((i, 1, 0.));
    }
    // Add edges from nodes in nodes_2 -> nodes_1 + nodes_2 + 1(end)
    for i in (nodes_1 + 1)..(nodes_1 + nodes_2 + 1) {
        graph[i].push((nodes_1 + nodes_2 + 1, 1, 0.));
    }
    let mut mcf = MinCostFlow::new(graph);
    // We have trivial solution res = 0, where there is no flow.
    // Thus, it is garanteed that the minimum cost flow exists.
    let res = mcf.min_flow_nolimit(start, end);
    debug!("maximum weight matching:{}", -res);
    assert!(res <= 0.);
    if res < 0. {
        mcf.enumerate_edge_used()
            .into_iter()
            .filter(|&(s, e)| s != 0 && e != nodes_1 + nodes_2 + 1)
            .map(|(s, e)| {
                if s < e {
                    (s - 1, e - nodes_1 - 1)
                } else {
                    (s - nodes_1 - 1, e - 1)
                }
            })
            .collect()
    } else {
        // No marging.
        vec![]
    }
}

#[derive(Debug, Clone)]
struct Edge {
    capacity: i64,
    cost: f64,
    to: usize,
    // If edge is u-> edge{to:v, rev:i},
    // then, the reverse edge would be edges[v][rev]
    rev: usize,
    // If this edge is the one in the original graph, true. otherwize false.
    is_original: bool,
}

impl Edge {
    fn new(to: usize, capacity: i64, cost: f64, rev: usize, b: bool) -> Self {
        Self {
            to,
            capacity,
            rev,
            cost,
            is_original: b,
        }
    }
}

const BIG: f64 = 100_000_000_000.;
#[derive(Debug, Clone)]
struct MinCostFlow {
    edges: Vec<Vec<Edge>>,
    size: usize,
}

impl MinCostFlow {
    // (index, capacity, cost)
    fn new(graph: Vec<Vec<(usize, i64, f64)>>) -> Self {
        let size = graph.len();
        let mut edges = vec![vec![]; size];
        for (from, targets) in graph.iter().enumerate() {
            for &(to, cap, cost) in targets {
                let rev_from = edges[to].len();
                let rev_to = edges[from].len();
                edges[from].push(Edge::new(to, cap, cost, rev_from, true));
                edges[to].push(Edge::new(from, 0, -cost, rev_to, false));
            }
        }
        Self { edges, size }
    }
    fn minimum_cost_flow(&mut self, start: usize, end: usize) -> Option<(u64, f64)> {
        // Bellman-Ford
        // The index of the edges used to get the score.
        // For example, if pred[to] = k,
        // then, the edge is edge = self.edges[to][k] and edge.to is the `from` node,
        // and the focal edges could get by self.edges[from][edge.rev].
        let mut pred: Vec<_> = vec![self.size + 1; self.size];
        let mut dist: Vec<_> = vec![BIG; self.size];
        dist[start] = 0.;
        for _ in 0..self.size - 1 {
            for (from, edges) in self.edges.iter().enumerate() {
                if dist[from] == BIG {
                    continue;
                }
                for edge in edges.iter().filter(|e| e.capacity > 0) {
                    let alternative_path = edge.cost + dist[from];
                    if alternative_path < dist[edge.to] {
                        dist[edge.to] = alternative_path;
                        pred[edge.to] = edge.rev;
                    }
                }
            }
        }
        if dist[end] == BIG {
            None
        } else {
            let mut cost = 0.;
            let mut temp = end;
            while temp != start {
                //eprint!("{}->", temp);
                let rev_edge_idx = pred[temp];
                let from_edge_idx = self.edges[temp][rev_edge_idx].rev;
                let next_node = self.edges[temp][rev_edge_idx].to;
                self.edges[temp][rev_edge_idx].capacity += 1;
                self.edges[next_node][from_edge_idx].capacity -= 1;
                cost += self.edges[next_node][from_edge_idx].cost;
                temp = next_node;
            }
            // eprintln!("cost:{}", cost);
            Some((1, cost))
        }
    }
    // Return minimum cost flow. There' no requirement on the total flow.
    // It is sometimes non-trivial because
    // there can be "negative cost flow."
    fn min_flow_nolimit(&mut self, start: usize, end: usize) -> f64 {
        let mut cost = 0.;
        while let Some((_f, c)) = self.minimum_cost_flow(start, end) {
            cost += c;
        }
        // The cost would be negative
        cost
    }
    // Enumerate all of used edges. It can be computed by
    fn enumerate_edge_used(&self) -> Vec<(usize, usize)> {
        self.edges
            .iter()
            .enumerate()
            .flat_map(|(idx, edges)| {
                edges
                    .iter()
                    .filter(|e| e.is_original)
                    .filter(|e| e.capacity == 0)
                    .map(|e| (idx, e.to))
                    .collect::<Vec<_>>()
            })
            .collect()
    }
}
