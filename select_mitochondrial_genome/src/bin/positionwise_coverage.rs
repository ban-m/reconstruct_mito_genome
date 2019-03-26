extern crate rayon;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::Result;
use std::path::Path;
#[derive(Debug)]
struct Coverage {
    coverage: Vec<u16>,
}

impl Coverage {
    fn convert_into_intervals(mappings: Vec<(usize, usize)>) -> Vec<(usize, i8)> {
        let mut intervals = Vec::with_capacity(mappings.len() * 2);
        for (start, end) in mappings {
            intervals.push((start, 1));
            intervals.push((end, -1));
        }
        intervals.sort_by_key(|e| e.0);
        intervals
    }
    fn new(coverage: Vec<u16>) -> Self {
        Coverage { coverage }
    }
    fn calc_positionwise_coverage(intervals: Vec<(usize, i8)>, length: usize) -> Vec<u16> {
        let (mut current_position, mut coverage) = (0, 0);
        let mut poswise = vec![];
        for &(pos, coverage_change) in &intervals {
            if pos != current_position {
                assert!(pos > current_position);
                (current_position..pos).for_each(|_| poswise.push(coverage));
                current_position = pos;
            }
            if coverage_change < 0 {
                coverage -= 1;
            } else {
                coverage += 1;
            }
        }
        (current_position..length).for_each(|_| poswise.push(coverage));
        poswise
    }
    fn from_interval(mappings: Vec<(usize, usize)>, length: usize) -> Self {
        let intervals = Self::convert_into_intervals(mappings);
        let coverages = Self::calc_positionwise_coverage(intervals, length);
        Coverage::new(coverages)
    }
    fn print(&self) {
        if !self.coverage.is_empty() {
            for i in 0..self.coverage.len() - 1 {
                print!(",{}", self.coverage[i]);
            }
            println!("{}", self.coverage.last().unwrap());
        }
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let paf = paf_open(&args[1])?;
    let coverage_of_each_read = positionwise_coverage(paf);
    for (id, coverage) in coverage_of_each_read {
        print!("{}", id);
        coverage.print();
    }
    Ok(())
}

fn positionwise_coverage(paf: String) -> HashMap<String, Coverage> {
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
        .map(|(id, mappings)| (id, Coverage::from_interval(mappings.0, mappings.1)))
        .collect()
}

// Return (read id of query, start position of query, end position of query, query length) and the
// same information of read target.
fn parse(line: &str) -> Option<(String, usize, usize, usize, String, usize, usize, usize)> {
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

fn paf_open(file: &str) -> Result<String> {
    let mut paf_raw: String = String::with_capacity(10_000_000_000);
    let mut paf = File::open(&Path::new(file))?;
    paf.read_to_string(&mut paf_raw)?;
    Ok(paf_raw)
}

#[test]
fn test_convert_into_intervals() {
    let mappings = vec![(0, 3), (1, 10), (2, 5)];
    let res = Coverage::convert_into_intervals(mappings);
    assert_eq!(
        res,
        vec![(0, 1), (1, 1), (2, 1), (3, -1), (5, -1), (10, -1)]
    );
    let mappings = vec![(0, 1), (2, 3), (4, 5)];
    let res = Coverage::convert_into_intervals(mappings);
    assert_eq!(res, vec![(0, 1), (1, -1), (2, 1), (3, -1), (4, 1), (5, -1)]);
    let mappings = vec![(0, 10), (1, 8), (3, 7), (5, 6)];
    let res = Coverage::convert_into_intervals(mappings);
    assert_eq!(
        res,
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
#[test]
fn test_summing_up() {
    let res = vec![(0, 1), (10, -1)];
    let cov = Coverage::calc_positionwise_coverage(res,10);
    assert_eq!(cov, vec![1;10]);
    let res = vec![(0, 1), (2, -1), (3, 1), (4, -1)];
    let cov = Coverage::calc_positionwise_coverage(res,4);
    assert_eq!(cov, vec![1,1,0,1]);
    let mappings = vec![(0, 3), (1, 10), (2, 5)];
    let res = Coverage::convert_into_intervals(mappings);
    let cov = Coverage::calc_positionwise_coverage(res,10);
    assert_eq!(cov, vec![1,2,3,2,2,1,1,1,1,1]);
    let mappings = vec![(0, 1), (2, 3), (4, 5)];
    let res = Coverage::convert_into_intervals(mappings);
    let cov = Coverage::calc_positionwise_coverage(res,5);
    assert_eq!(cov, vec![1,0,1,0,1]);
    let mappings = vec![(0, 10), (1, 8), (3, 7), (5, 6)];
    let res = Coverage::convert_into_intervals(mappings);
    let cov = Coverage::calc_positionwise_coverage(res,10);
    assert_eq!(cov,vec![1,2,2,3,3,4,3,2,1,1]);
}
