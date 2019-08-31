extern crate bio;
extern crate handmade_bloom_filter;
extern crate rayon;
use bio::io;
use handmade_bloom_filter::UpperBoundBFFactory;
use handmade_bloom_filter::HUGE_MODULO;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let t: usize = args[2].parse().unwrap();
    let k: usize = args[3].parse().unwrap();
    let assignment: HashMap<_, _> = open_assignment(&args[4])?;
    let input: Vec<_> = io::fastq::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .map(|e| e.seq().to_vec())
        .collect();
    let bf = UpperBoundBFFactory::default()
        .k(k)
        .number_of_hash(7)
        .modulo(HUGE_MODULO)
        .add_dataset(&input)
        .finalize_par(t);
    eprintln!("{}", bf);
    let input: Vec<_> = io::fastq::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| assignment.contains_key(e.id()))
        .collect();
    use std::io::{BufWriter, Write};
    let mut wtr = BufWriter::new(std::io::stdout());
    for record in input {
        let line = bf.upper_bound_at_each_position(record.seq());
        let id = record.id();
        writeln!(&mut wtr, "{},{}", id, to_string(&line))?;
    }
    Ok(())
}

fn mean_sd(xs: &Vec<u16>) -> (u64, u64) {
    let mut sum = 0;
    let mut sumsq = 0;
    let len = xs.len() as u64;
    for x in xs {
        sum += *x as u64;
        sumsq += (*x as u64).pow(2);
    }
    let ave = sum / len;
    let sd = ((sumsq / len - ave * ave) as f64).sqrt().floor() as u64;
    (ave, sd)
}

fn to_string(cov: &[u16]) -> String {
    let cov: Vec<_> = cov.iter().map(|e| format!("{}", e)).collect();
    cov.join(",")
}

fn open_assignment(file: &str) -> std::io::Result<HashMap<String, String>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;
    let mut res = HashMap::new();
    for line in BufReader::new(File::open(&Path::new(file))?)
        .lines()
        .filter_map(|e| e.ok())
    {
        let mut contents = line.split('\t');
        let id = contents.next().unwrap().to_string();
        let chrtype = contents.next().unwrap();
        let chrtype = match chrtype.parse::<u8>() {
            Ok(_) => "genomic".to_string(),
            Err(_) => chrtype.to_string(),
        };
        res.insert(id, chrtype);
    }
    Ok(res)
}
