extern crate bio_utils;
extern crate rusty_sandbox;
use std::collections::HashMap;

fn calc_connected_component(graph: &[Vec<usize>]) {
    // If more than 1, the node has been already arrived.
    // A node without any edges is consided 0-component, as it does not exists, in fact.
    let mut nodes: Vec<usize> = graph
        .iter()
        .map(|e| if e.is_empty() { 0 } else { 1 })
        .collect();
    let mut current_color = 2;
    let mut idx = 0;
    while idx < nodes.len() {
        eprintln!("Searching {}-th component... Index:{}", current_color, idx);
        // Depth first search.
        let mut stack = vec![];
        stack.push(idx);
        while !stack.is_empty() {
            let child = stack.pop().unwrap();
            if nodes[child] != current_color {
                nodes[child] = current_color;
                stack.extend(graph[child].iter());
            }
        }
        eprintln!("Finishing searchig this color....");
        // Find next starting point.
        while idx < nodes.len() && nodes[idx] != 1 {
            idx += 1;
        }
        current_color += 1;
    }
    eprintln!("Finishing coloring graph");
    let mut summary = HashMap::new();
    for node in 0..nodes.len() {
        // (# nodes, # edges)
        let entry = summary.entry(nodes[node]).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += graph[node].len();
    }
    println!("Component\tNodes\tEdges");
    for component in 2..current_color {
        let (nodes, edges) = summary[&component];
        println!("{}\t{}\t{}", component, nodes, edges);
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (graph, _id_to_index) = rusty_sandbox::construct_graph(&args[1])?;
    let mut directed_graph = vec![vec![]; graph.len()];
    for (from, edges) in graph.into_iter().enumerate() {
        for (to, _, _) in edges {
            directed_graph[from].push(to);
            directed_graph[to].push(from);
        }
    }
    eprintln!("Constructed Graph");
    calc_connected_component(&directed_graph);
    Ok(())
}
