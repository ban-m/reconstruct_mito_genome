use std::path::Path;
use std::io::{BufReader,BufRead};
use std::collections::HashSet;
use std::fs::File;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let first_id:HashSet<_> = BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e|e.ok())
        .collect();
    let second_id:HashSet<_> = BufReader::new(File::open(&Path::new(&args[2]))?)
        .lines()
        .filter_map(|e|e.ok())
        .collect();
    println!("ReadID\tInFirst\tIsInSecond");
    for id in &first_id{
        if second_id.contains(id){
            println!("{}\t{}\t{}",id,1,1);
        }else{
            println!("{}\t{}\t{}",id,1,0);
        }
    }
    for id in &second_id{
        if !first_id.contains(id){
            println!("{}\t{}\t{}",id,0,1);
        }
    }
    Ok(())
}
