use std::io::Result;
use std::io::{BufReader,BufRead};
use std::fs::File;
fn main()-> Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let summary = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e|e.ok())
        .map(|e| {
            let contents:Vec<_> = e.split(',').collect();
            let data:Vec<u32> = contents[1..].iter()
                .filter_map(|e|e.parse().ok()).collect();
            (contents[0].to_owned(),data)
        })
        .map(|(id,data)| (id,summarize(data)));
    for (id, (cov, p)) in summary {
            println!("{}\t{}\t{}",id,cov,p);
    }
    Ok(())
}

fn summarize(data:Vec<u32>)->(f64,f64){
    let sum:u32 = data.iter().sum();
    let _mean = sum as f64 / data.len() as f64;
    (0.,0.)
}
