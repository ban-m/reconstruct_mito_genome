extern crate bio_utils;
extern crate rusty_sandbox;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (graph, _id_to_index) = rusty_sandbox::construct_graph(&args[1])?;
    let graph:Vec<Vec<usize>> = graph.into_iter().map(|xs| xs.into_iter().map(|e|e.0).collect()).collect();
    eprintln!("Constructed Graph");
    let components = rusty_sandbox::calc_connected_component(&graph);
    println!("Component\tNodes\tEdges");
    let mut count =0;
    for component in components{
        let nodes = component.len();
        let edges:usize = component
            .iter()
            .map(|&node| graph[node].iter().filter(|e| component.contains(e)).count())
            .sum();
        println!("{}\t{}\t{}",count, nodes, edges);
        count += 1;
    }
    Ok(())
}
