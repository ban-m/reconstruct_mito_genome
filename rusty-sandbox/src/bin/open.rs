extern crate rayon;
extern crate rusty_sandbox;
use rayon::prelude::*;
use std::time;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let now = time::Instant::now();
    let result = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    // let result = open(&args[1])?;
    let elapsed = time::Instant::now() - now;
    println!("Open:{:?}", elapsed);
    let now = time::Instant::now();
    let num: usize = result
        .into_par_iter()
        .flat_map(|e| e.1)
        .map(|e| e.cigar().len())
        .sum();
    let elapsed = time::Instant::now() - now;
    println!("Collect:{:?}", elapsed);
    println!("{}", num);
    Ok(())
}
