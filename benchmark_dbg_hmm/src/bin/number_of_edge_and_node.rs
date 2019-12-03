extern crate dbg_hmm;
extern crate edlib_sys;
extern crate env_logger;
extern crate rand;
extern crate log;
use dbg_hmm::gen_sample::*;
use dbg_hmm::*;
use rand::{rngs::StdRng, SeedableRng};
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 150;
    let num_seq = 250;
    let mut rng: StdRng = SeedableRng::seed_from_u64(12_121_899_892);
    let template: Vec<_> = generate_seq(&mut rng, len);
    let data: Vec<Vec<_>> = (0..num_seq)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let k = 6;
    let mut f = Factory::new();
    println!("Coverage\tNodes\tEdges\tWeight\tMethod");
    for i in 1..num_seq {
        let ms: Vec<_> = data[..i].iter().map(|e| e.as_slice()).collect();
        let w = vec![1.; i];
        let m = f.generate_with_weight(&ms, &w, k);
        println!("{}\t{}\t{}\t{}\tW", i, m.node_num(), m.edge_num(), m.weight());
        let m = f.generate_from_ref(&ms, k);
        println!("{}\t{}\t{}\t{}\tH", i, m.node_num(), m.edge_num(), m.weight());
    }
}
