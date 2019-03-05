extern crate bio;
use std::io::{BufRead,BufReader};
use std::fs::File;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let coverage = get_coverage(&args[1])?;
    for (id, cov) in coverage{
        println!("{}\t{}",id,cov);
    }
    Ok(())
}
use std::collections::HashMap;

fn get_coverage(path:&str)->std::io::Result<Vec<(String,u32)>>{
    let mut seq:HashMap<_,_> = HashMap::new();
    for (id1,id2) in BufReader::new(File::open(&std::path::Path::new(path))?)
        .lines()
        .filter_map(|e|e.ok())
        .map(|e|{
            let e:Vec<_> = e.split('\t').collect();
            (e[0].to_string(), e[5].to_string())}){
            let cov = seq.entry(id1).or_insert(0);
            *cov +=1;
            let cov = seq.entry(id2).or_insert(0);
            *cov +=1;
        }
    Ok(seq.into_iter().collect())
}

