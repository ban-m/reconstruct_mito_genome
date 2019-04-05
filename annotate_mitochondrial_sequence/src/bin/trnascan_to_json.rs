extern crate rayon;
extern crate bio;
extern crate annotate_mitochondrial_sequence;
use std::io::{BufReader,BufRead};
use std::fs::File;
use std::path::Path;
use annotate_mitochondrial_sequence::Map;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let input:Vec<_> = BufReader::new(File::open(&Path::new(&args[1]))?)
        .lines()
        .filter_map(|e|e.ok())
        .skip_while(|e|!e.starts_with("----"))
        .skip(1)// delimiter itself.
        .map(|e|TRNAScan::new(&e))
        .map(|e| e.to_map())
        .collect();
    println!("{}",serde_json::ser::to_string_pretty(&input).unwrap());
    Ok(())
}

#[derive(Debug)]
struct TRNAScan{
    name:String,
    start:usize,
    end:usize,
    strand:i8,
    protein:String,
    score:f64,
}

impl TRNAScan{
    fn new(line:&str)-> Self{
        let e:Vec<_> = line.split('\t').collect();
        let name = e[0].trim().to_string();
        let pos1 = match e[2].trim().parse(){
            Ok(res)=>res,
            Err(why) => {eprintln!("2 -> {}, {}, {:?}",line,e[2],why);
                         panic!();},
        };
        let pos2 = match e[3].trim().parse(){
            Ok(res)=>res,
            Err(why) => {eprintln!("3 -> {}, {},{:?}",line,e[3],why);
                         panic!();},
        };

        let (start,end,strand) = Self::determine(pos1,pos2);
        let protein = e[4].to_string();
        let score = match e[8].trim().parse::<f64>(){
            Ok(res)=>res / 100.,
            Err(why) => {eprintln!("8 -> {}, {}, {:?}",line,e[8],why);
                         panic!();},
        };
        TRNAScan{
            name,
            start,
            end,
            strand,
            protein,
            score
        }
    }
    fn determine(pos1:usize,pos2:usize)->(usize,usize,i8){
        if pos1 < pos2 {
            (pos1,pos2,1)
        }else{
            (pos2,pos1,-1)
        }
    }
    fn to_map(&self)->Map{
        Map::new(&self.name,
                 true,
                 self.strand == 1,
                 self.start as u64,
                 self.end as u64,
                 "tRNA-SE",
                 &self.protein,
                 self.score)
    }
}
