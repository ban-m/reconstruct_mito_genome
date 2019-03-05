extern crate rand;
use rand::{thread_rng};
use rand::prelude::*;
use std::collections::HashSet;
pub mod findunion;
pub fn pick_first_cycle(nodes_size:usize,edges:&[(usize,usize)])->(usize,Vec<usize>){
    let mut fu = findunion::FindUnion::new(nodes_size);
    let mut order:Vec<_> = (0..edges.len()).collect();
    let mut rng = thread_rng();
    let mut pushed_nodes = HashSet::with_capacity(10000);
    let mut picked_edges = 0;
    order.shuffle(&mut rng);
    for edge_idx in order{
        picked_edges += 1;
        let (u,v) = edges[edge_idx];
        pushed_nodes.insert(u);
        pushed_nodes.insert(v);
        if fu.same(u,v).unwrap(){
            // a cycle would be created by adding (u,v) in subgraph.
            let res = pushed_nodes.into_iter().filter(|&e|fu.same(e,u).unwrap_or(false))
                .collect();
            return (picked_edges,res)
        }else{
            eprintln!("Unite {} and {}",u,v);
            fu.unite(u,v).unwrap();
            eprintln!("{:?}",fu.same(u,v));
        }
    }
    (picked_edges,vec![])
}

#[test]
fn pick_first_cycle_test(){
    let edges:Vec<_> = (0..10).map(|e|(e,(e+1) % 10)).collect();
    let (picknum,res) = pick_first_cycle(11,&edges);
    eprintln!("{:?},{}",res,picknum);
    assert_eq!(picknum,10);
}
