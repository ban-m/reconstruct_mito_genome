use std::fs::File;
use std::io::{Read,BufRead,BufReader};
use std::collections::HashMap;
use std::path::Path;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let input:String = open_file(&args[1])?;
    let assignment:HashMap<_,_> = open_assignment(&args[2])?;
    let (mut mito,mut chloro, mut genome, mut genome_low, mut unmap) = (0,0,0,0,0);
    for line in input.lines(){
        let mut content = line.split(',');
        let id = match content.next(){
            Some(res) => res,
            None => continue,
        };
        if assignment[id] == "genomic"{
            let average = average(content);
            if genome < 5 && average > 500 {
                genome += 1;
                println!("{}",line);
            }else if genome_low < 5 && average < 500 {
                genome_low += 1;
                println!("{}",line);
            }
        }
        if assignment[id] == "mitochondria" && mito < 5 {
            mito += 1;
            println!("{}",line);
        }else if assignment[id] == "chloroplast" && chloro < 5{
            chloro += 1;
            println!("{}",line);
        }else if assignment[id] == "*" && unmap < 5 {
            unmap += 1;
            println!("{}",line);
        }
        if mito > 5 && chloro > 5 && genome > 5 && genome_low > 5 && unmap > 5 {
            break;
        }
    }
    Ok(())
}

fn open_assignment(file:&str)->std::io::Result<HashMap<String,String>>{
    let mut res = HashMap::new();
    for line in BufReader::new(File::open(&Path::new(file))?)
        .lines()
        .filter_map(|e|e.ok())
    {
        let mut contents = line.split('\t');
        let id = contents.next().unwrap().to_string();
        let chrtype = contents.next().unwrap();
        let chrtype = match chrtype.parse::<u8>(){
            Ok(_) => "genomic".to_string(),
            Err(_) => chrtype.to_string()
        };
        res.insert(id,chrtype);
    }
    Ok(res)
}

fn open_file(file:&str)->std::io::Result<String>{
    let mut input = String::new();
    let mut file = File::open(&Path::new(file))?;
    file.read_to_string(&mut input)?;
    Ok(input)
}

fn average<'a>(content: std::str::Split<'a,char>)->u64{
    let (sum,len) = content.filter_map(|e|e.parse::<u64>().ok())
        .fold((0,0),|(acc,len),x|(acc + x , len + 1));
    sum / len
}
