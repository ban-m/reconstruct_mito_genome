use std::collections::HashSet;
use std::io::{BufReader,BufRead};
use std::fs::File;
use std::path::Path;
fn main()-> std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let read_ids_1 = get_read_ids(&args[1])?;
    let read_ids_2 = get_read_ids(&args[2])?;
    let first_size = read_ids_1.len();
    let second_size = read_ids_2.len();
    let intersect_size = read_ids_1.intersection(&read_ids_2).count();
    println!("first bulk\t{}",read_ids_1.len());
    println!("second bulk\t{}",read_ids_2.len());
    println!("first only\t{}",first_size - intersect_size);
    println!("second only\t{}",second_size - intersect_size);
    println!("intersect_size\t{}",intersect_size);
    Ok(())
}

fn get_read_ids(file:&str)->std::io::Result<HashSet<String>>{
    Ok(
        BufReader::new(File::open(&Path::new(file))?).lines()
            .filter_map(|e|e.ok())
            .filter(|e|!e.is_empty())
            .collect())
}
