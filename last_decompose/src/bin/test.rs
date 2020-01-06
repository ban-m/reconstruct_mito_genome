extern crate last_decompose;
fn main() {
    let nodes1 = 3;
    let nodes2 = 3;
    let graph = vec![
        vec![(0, 3.), (1, 5.), (2, 1.)],
        vec![(0, 7.), (1, 8.), (2, 6.)],
        vec![(0, 1.), (1, 1.), (2, 5.)],
    ];
    let res = last_decompose::bipartite_matching::maximum_weight_matching(nodes1, nodes2, &graph);
    eprintln!("{:?}", res);
}
