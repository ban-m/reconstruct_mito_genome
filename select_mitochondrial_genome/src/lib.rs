extern crate rayon;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::Result;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn paf_open(file: &str) -> Result<String> {
    let mut paf_raw: String = String::with_capacity(10_000_000_000);
    let mut paf = File::open(&Path::new(file))?;
    paf.read_to_string(&mut paf_raw)?;
    Ok(paf_raw)
}

// Return (read id of query, start position of query, end position of query, query length) and the
// same information of read target.
pub fn parse(line: &str) -> Option<(String, usize, usize, usize, String, usize, usize, usize)> {
    let contents: Vec<&str> = line.split('\t').collect();
    let query_length: usize = contents[1].parse().ok()?;
    let query_start: usize = contents[2].parse().ok()?;
    let query_end: usize = contents[3].parse().ok()?;
    let target_length: usize = contents[6].parse().ok()?;
    let target_start: usize = contents[7].parse().ok()?;
    let target_end: usize = contents[8].parse().ok()?;
    let _alignment_block_length: usize = contents[10].parse().ok()?;
    Some((
        contents[0].to_string(),
        query_start,
        query_end,
        query_length,
        contents[5].to_string(),
        target_start,
        target_end,
        target_length,
    ))
}

#[derive(Debug)]
pub struct Interval {
    inner: Vec<(usize, i8)>,
    length: usize,
}

impl Interval {
    pub fn new(mappings: &[(usize, usize)], length: usize) -> Self {
        let mut intervals = Vec::with_capacity(mappings.len() * 2);
        for &(start, end) in mappings {
            intervals.push((start, 1));
            intervals.push((end, -1));
        }
        intervals.sort_by_key(|e| e.0);
        Interval {
            inner: intervals,
            length: length,
        }
    }
    pub fn from_raw(inner: &[(usize, i8)], length: usize) -> Self {
        Interval {
            inner: inner.to_vec(),
            length: length,
        }
    }
    pub fn length(&self) -> usize {
        self.length
    }
    pub fn inner(&self) -> &Vec<(usize, i8)> {
        &self.inner
    }
}

pub fn paf_to_intervals(paf: String) -> Vec<(String, Interval)> {
    let mut summary_of_each_read: HashMap<String, (Vec<(usize, usize)>, usize)> = HashMap::new();
    for line in paf.lines() {
        if let Some((read1_id, read1_s, read1_e, length1, read2_id, read2_s, read2_e, length2)) =
            parse(line)
        {
            let entry = summary_of_each_read
                .entry(read1_id)
                .or_insert((vec![], length1));
            (entry.0).push((read1_s, read1_e));
            let entry = summary_of_each_read
                .entry(read2_id)
                .or_insert((vec![], length2));
            (entry.0).push((read2_s, read2_e));
        }
    }
    summary_of_each_read
        .into_par_iter()
        .map(|(id, mappings)| (id, Interval::new(&mappings.0, mappings.1)))
        .collect()
}

use std::collections::HashSet;
fn get_ids(file: &str) -> HashSet<String> {
    use std::io::{BufRead, BufReader};
    BufReader::new(std::fs::File::open(&std::path::Path::new(file)).unwrap())
        .lines()
        .skip(1)
        .filter_map(|e| e.ok())
        .collect()
}

pub fn paf_file_to_intervals_with_id(paf: &str, ids: &str) -> Vec<(String, Interval)> {
    eprintln!("Opening file");
    let ids = get_ids(ids);
    eprintln!("Opened id file");
    let mut summary_of_each_read: HashMap<String, (Vec<(usize, usize)>, usize)> = HashMap::new();
    let reader = BufReader::new(File::open(&Path::new(paf)).unwrap());
    // let mut line = String::new();
    // while reader.read_line(&mut line).unwrap() > 0 {
    //     if line.is_empty(){
    //         break;
    //     }
    for line in reader.lines().filter_map(|e|e.ok()){
        if let Some((read1_id, read1_s, read1_e, length1, read2_id, read2_s, read2_e, length2)) =
            parse(&line)
        {
            if !ids.contains(&read1_id) || !ids.contains(&read2_id){
                continue;
            }
            let entry = summary_of_each_read
                .entry(read1_id)
                .or_insert((vec![], length1));
            (entry.0).push((read1_s, read1_e));
            let entry = summary_of_each_read
                .entry(read2_id)
                .or_insert((vec![], length2));
            (entry.0).push((read2_s, read2_e));
        }
    }
    eprintln!("Finish read file. {} reads.",summary_of_each_read.len());
    summary_of_each_read
        .into_par_iter()
        .map(|(id, mappings)| (id, Interval::new(&mappings.0, mappings.1)))
        .collect()
}

pub fn paf_file_to_intervals(paf: &str) -> Vec<(String, Interval)> {
    let mut summary_of_each_read: HashMap<String, (Vec<(usize, usize)>, usize)> = HashMap::new();
    let mut reader = BufReader::new(File::open(&Path::new(paf)).unwrap());
    let mut line = String::new();
    while reader.read_line(&mut line).unwrap() > 0 {
        if let Some((read1_id, read1_s, read1_e, length1, read2_id, read2_s, read2_e, length2)) =
            parse(&line)
        {
            let entry = summary_of_each_read
                .entry(read1_id)
                .or_insert((vec![], length1));
            (entry.0).push((read1_s, read1_e));
            let entry = summary_of_each_read
                .entry(read2_id)
                .or_insert((vec![], length2));
            (entry.0).push((read2_s, read2_e));
        }
    }
    summary_of_each_read
        .into_par_iter()
        .map(|(id, mappings)| (id, Interval::new(&mappings.0, mappings.1)))
        .collect()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    #[test]
    fn test_convert_into_intervals() {
        let mappings = vec![(0, 3), (1, 10), (2, 5)];
        let res = Interval::new(&mappings, 10);
        assert_eq!(
            res.inner,
            vec![(0, 1), (1, 1), (2, 1), (3, -1), (5, -1), (10, -1)]
        );
        let mappings = vec![(0, 1), (2, 3), (4, 5)];
        let res = Interval::new(&mappings, 5);
        assert_eq!(
            res.inner,
            vec![(0, 1), (1, -1), (2, 1), (3, -1), (4, 1), (5, -1)]
        );
        let mappings = vec![(0, 10), (1, 8), (3, 7), (5, 6)];
        let res = Interval::new(&mappings, 10);
        assert_eq!(
            res.inner,
            vec![
                (0, 1),
                (1, 1),
                (3, 1),
                (5, 1),
                (6, -1),
                (7, -1),
                (8, -1),
                (10, -1)
            ]
        );
    }
}
