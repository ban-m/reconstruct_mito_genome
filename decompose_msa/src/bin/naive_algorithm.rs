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
    for (i, e) in transpose
        .into_iter()
        .map(|column| entropy(&column))
        .enumerate()
    {
        println!("{},{}", i, e);
    }
    Ok(())
}

fn entropy(column: &[u8]) -> f64 {
    let count = column.iter().fold([0; 4], |mut res, base| {
        match base {
            b'A' => res[0] += 1,
            b'C' => res[1] += 1,
            b'G' => res[2] += 1,
            b'T' => res[3] += 1,
            _ => {}
        };
        res
    });
    let sum: i32 = count.iter().sum();
    count
        .iter()
        .map(|&count| {
            if count == 0 {
                0.
            } else {
                let f = count as f64 / sum as f64;
                -f * f.ln()
            }
        })
        .sum()
}
