extern crate rand;
use rand::{seq::*, thread_rng};
fn main() {
    let mut rng = thread_rng();
    let seq = b"123456";
    println!(
        "{}",
        index::sample(&mut rng, seq.len(), seq.len())
            .into_iter()
            .map(|idx| seq[idx] as char)
            .collect::<String>()
    );
}
