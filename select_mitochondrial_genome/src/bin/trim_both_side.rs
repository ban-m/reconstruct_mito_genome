use std::fs::File;
use std::path::Path;
use std::io::Read;
use std::io::{BufWriter};
fn main() -> std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let input = open_file(&args[1])?;
    let _wtr = BufWriter::new(std::io::stdout());
    for line in input.lines(){
        let mut contents = line.split(',');
        let _id = contents.next().unwrap();
        let coverage:Vec<u32> = contents.filter_map(|e| e.parse().ok()).collect();
        let len = coverage.len();
        let _coverage:Vec<u32> = coverage.into_iter().skip(len*5/100).take(len*95/100).collect();
    }
    Ok(())
    
}

fn open_file(file:&str)->std::io::Result<String>{
    let mut file = File::open(&Path::new(file))?;
    let mut input = String::new();
    file.read_to_string(&mut input)?;
    Ok(input)
}
