use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let answer: HashMap<String, usize> = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let e: Vec<_> = e.split('\t').collect();
            let cl = e[0].parse().ok()?;
            let id = e[1].to_string();
            Some((id, cl))
        })
        .collect();
    let max = *answer.values().max().unwrap();
    let mut result = vec![vec![0; max as usize]; max as usize];
    let mut count = vec![0; max as usize];
    for line in BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
    {
        let e: Vec<_> = line.split('\t').collect();
        let cl = match e[0].parse::<usize>().ok() {
            Some(res) => res,
            None => continue,
        };
        let answer = match answer.get(e[1]) {
            Some(res) => *res,
            None => continue,
        };
        result[answer][cl] += 1;
        count[answer] += 1;
    }
    for (res, count) in result.into_iter().zip(count) {
        let res: Vec<String> = res.iter().map(|e| format!("{}", e)).collect();
        let res = res.join("\t");
        eprintln!("{}\t{}", res, count);
    }
    Ok(())
}