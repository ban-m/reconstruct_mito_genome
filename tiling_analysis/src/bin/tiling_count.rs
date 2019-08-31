use std::collections::HashMap;
use std::io::Read;

fn count_line(line:&str)-> Vec<usize>{
    let encoded = match line.split(':').nth(1){
        Some(res) => res,
        None => return vec![],
    };
    let mut units:Vec<_> = encoded
        .split_whitespace()
        .filter_map(|e|e.split('-').nth(0).and_then(|e|e.parse::<usize>().ok()))
        .collect();
    units.sort();
    units.dedup();
    units
}

fn count(file:&str)->std::io::Result<Vec<usize>>{
    let mut file = std::fs::File::open(file)?;
    let mut input = String::new();
    file.read_to_string(&mut input)?;
    Ok(input.lines().flat_map(count_line).collect())
}

fn count_from_file(file:&str)->std::io::Result<HashMap<usize, usize>>{
    let mut result = HashMap::new();
    for unit in count(file)?{
        let count = result.entry(unit).or_insert(0);
        *count += 1;
    }
    Ok(result)
}

fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    for (unit, count) in count_from_file(&args[1])?{
        println!("{}\t{}", unit,count);
    }
    Ok(())
}
