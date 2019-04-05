use std::fs::File;
use std::io::{Read,BufRead,BufReader,BufWriter,Write};
use std::collections::HashMap;
use std::path::Path;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let starts = open_file(&args[1])?;
    let stops = open_file(&args[2])?;
    let mut starts_out = BufWriter::new(File::create(&Path::new(&args[4]))?);
    let mut stops_out = BufWriter::new(File::create(&Path::new(&args[5]))?);
    let assignment:HashMap<_,_> = open_assignment(&args[3])?;
    let (mut mito,mut chloro, mut genome, mut unmap) = (0,0,0,0);
    for (start,stop) in starts.lines().zip(stops.lines()){
        let (mut start_con,mut stop_con) = (start.split(','),stop.split(','));
        let (id1,id2) = match (start_con.next(), stop_con.next()){
            (Some(re1),Some(re2)) => (re1,re2),
            _ => continue,
        };
        if !is_ok_data(start_con,stop_con){
            continue;
        }
        assert_eq!(id1,id2);
        if assignment[id1] == "genomic" && genome < 5 {
            genome += 1;
            writeln!(&mut starts_out, "{}",start)?;
            writeln!(&mut stops_out, "{}",stop)?;
        }else if assignment[id1] == "mitochondria" && mito < 5 {
            mito += 1;
            writeln!(&mut starts_out, "{}",start)?;
            writeln!(&mut stops_out, "{}",stop)?;
        }else if assignment[id1] == "chloroplast" && chloro < 5{
            chloro += 1;
            writeln!(&mut starts_out, "{}",start)?;
            writeln!(&mut stops_out, "{}",stop)?;
        }else if assignment[id1] == "*" && unmap < 5 {
            unmap += 1;
            writeln!(&mut starts_out, "{}",start)?;
            writeln!(&mut stops_out, "{}",stop)?;
        }
        if mito > 5 && chloro > 5 && genome > 5 && unmap > 5 {
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

fn is_ok_data<'a>(mut start: std::str::Split<'a,char>,mut stop: std::str::Split<'a,char>)->bool{
    // header is alredy trimed.
    let len:usize = match start.next().and_then(|e|e.parse().ok()){
        Some(res) => res,
        _ => return false
    };
    let len2:usize = match stop.next().and_then(|e|e.parse().ok()){
        Some(res) => res,
        _ => return false
    };
    len > 10 && len2 > 10
}
