use std::fs::File;
use std::io::{Read,Write};
use std::path::Path;
const LEN:usize = 1000;
// Used for t-SNE
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let input:String = open_file(&args[1])?;
    let mut output = vec![];
    for line in input.lines(){
        let (id,coverage) = parse(line);
        let normalized = normalize(&coverage);
        write!(&mut output,"{}",id)?;
        for cov in normalized{
            write!(&mut output,",{}",cov)?;
        }
        writeln!(&mut output,"")?;
    }
    println!("{}",String::from_utf8(output).unwrap());
    Ok(())
}

fn open_file(file:&str)->std::io::Result<String>{
    let mut input = String::new();
    let mut file = File::open(&Path::new(file))?;
    file.read_to_string(&mut input)?;
    Ok(input)
}

fn parse<'a>(line:&'a str)->(&'a str,Vec<u32>){
    let mut content = line.split(',');
    let id = content.next().unwrap();
    let cov:Vec<u32> = content.filter_map(|e|e.parse().ok()).collect();
    (id,cov)
}
// Rather smoothing.
fn normalize(cov:&[u32])->Vec<u32>{
    let mut result = vec![0;LEN];
    let len = cov.len();
    for i in 0..LEN{
        // calculate i-th bin;
        let mut length_of_ith_bin = 0;
        let mut sum_of_ith_bin = 0;
        let (start,end) = ((i*len)/LEN,((i+1)*len)/LEN);
        if start == end {
            for j in start..end+1{
                sum_of_ith_bin += cov[j];
                length_of_ith_bin +=1;
            }
        }else{
            for j in start..end{
                sum_of_ith_bin += cov[j];
                length_of_ith_bin +=1;
            }
        }
        result[i] = sum_of_ith_bin/length_of_ith_bin;
    }
    result
}
