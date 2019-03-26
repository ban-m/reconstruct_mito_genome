use std::path::Path;
use std::io::Read;
use std::fs::File;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let input = open_paf_file(&args[1])?;
    let graph = summarize(&input);
    println!("ID\tInDegree\tOutDegree");
    for (id,indegree,outdegree) in &graph{
        println!("{}\t{}\t{}",id,indegree,outdegree);
    }
    // let (mut max_in, mut sum_in, mut sumsq_in) = (0,0,0);
    // let (mut max_out, mut sum_out, sumsq_out) = (0,0,0);
   Ok(())
}

fn open_paf_file(file:&str)->std::io::Result<String>{
    let mut reader = File::open(&Path::new(file)).unwrap();
    let mut input = String::new();
    reader.read_to_string(&mut input)?;
    Ok(input)
}

fn summarize(_input:&str)->Vec<(String,usize,usize)>{
    vec![]
}
