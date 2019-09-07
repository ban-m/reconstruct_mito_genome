use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let msa: Vec<Vec<_>> = BufReader::new(std::fs::File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .map(|e| e.as_bytes().to_vec())
        .collect();
    let transpose = {
        let mut res = vec![vec![0; msa.len()]; msa[0].len()];
        for i in 0..msa.len() {
            for j in 0..msa[i].len() {
                res[j][i] = msa[i][j];
            }
        }
        res
    };
    for i in 0..transpose.len(){
        for j in i..transpose.len(){
            println!("{},{},{}",hypergeom(&transpose[i],&transpose[j]));
        }
    }
    Ok(())
}

fn hypergeom(a:&[u8],b:&[u8])->f64{
    let large_n = a.len();
    let (large_k, major)
        = {
        let (a,t,g,c) = 
    };
}
