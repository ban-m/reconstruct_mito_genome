extern crate dbg_hmm;
extern crate env_logger;
use dbg_hmm::Factory;
use std::fs::File;
use std::io::{BufRead, BufReader};
fn main() {
    env_logger::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let lines: Vec<_> = BufReader::new(File::open(&args[1]).unwrap())
        .lines()
        .filter_map(|e| e.ok())
        .collect();
    let (data, ws): (Vec<Vec<u8>>, Vec<f64>) = lines
        .into_iter()
        .map(|e| {
            let mut e = e.split('\t');
            let seq: Vec<_> = e.next().unwrap().as_bytes().to_vec();
            let w: f64 = e.next().unwrap().parse().unwrap();
            //println!("{}\t{:.3}", String::from_utf8_lossy(&seq), w);
            (seq, w)
        })
        .unzip();
    let data: Vec<_> = data.iter().map(|e| e.as_slice()).collect();
    let mut f = Factory::new();
    let m = f.generate_with_weight(&data, &ws, 6, &mut vec![]);
    println!("{}", m);
}
