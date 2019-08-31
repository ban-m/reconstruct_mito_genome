extern crate rayon;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::Result;
use std::io::{BufRead, BufReader};
use std::path::Path;

const THRESHOLD: Option<usize> = None;
#[derive(Debug)]
pub struct ReadSummary {
    pub mean: f64,
    pub sd: f64,
    pub length: usize,
}

impl ReadSummary {
    pub fn new(sum: usize, sumsq: usize, length: usize) -> Self {
        let mean = sum as f64 / length as f64;
        let variance = sumsq as f64 / length as f64 - mean * mean;
        ReadSummary {
            mean: mean,
            sd: variance.sqrt(),
            length: length,
        }
    }
    pub fn summing_up(intervals: &Interval) -> (usize, usize) {
        let (mut current_position, mut coverage) = (0, 0);
        let (mut sum, mut sumsq) = (0, 0);
        for &(pos, coverage_change) in intervals.inner().iter() {
            if pos > current_position {
                (current_position..pos).for_each(|_| {
                    sum += coverage;
                    sumsq += coverage * coverage;
                });
                current_position = pos;
            }
            coverage = Self::update(coverage, coverage_change);
        }
        (sum, sumsq)
    }
    #[inline]
    pub fn update(current: usize, change: i8) -> usize {
        if change < 0 {
            current - 1
        } else {
            current + 1
        }
    }
    // Input:intervals Output: sum, sum of sqare, number of considered position.
    pub fn summing_up_trim(interval: &Interval, is_trim_on: usize) -> (usize, usize, usize) {
        // let threshold = Self::get_threshold(interval, is_trim_on);
        let (mut current_position, mut coverage) = (0, 0);
        let start = interval.length() * 3 / 100;
        let end = interval.length() * 97 / 100;
        let (mut sum, mut sumsq) = (0, 0);
        for &(pos, coverage_change) in interval
            .inner()
            .iter()
            .skip_while(|&(pos, _)| pos < &start)
            .take_while(|&(pos, _)| pos < &end)
        {
            if pos > current_position {
                (current_position..pos).for_each(|_| {
                    sum += coverage;
                    sumsq += coverage * coverage;
                });
                current_position = pos;
            }
            coverage = Self::update(coverage, coverage_change);
        }
        let len = end - start;
        eprintln!("{}->{} ({}% trimed)", interval.length(), len, is_trim_on);
        (sum, sumsq, len)
    }
    pub fn from_interval(interval: Interval) -> Self {
        if let Some(is_trim_on) = THRESHOLD {
            let (sum, sumsq, length) = Self::summing_up_trim(&interval, is_trim_on);
            ReadSummary::new(sum, sumsq, length)
        } else {
            let (sum, sumsq) = Self::summing_up(&interval);
            ReadSummary::new(sum, sumsq, interval.length())
        }
    }
}

pub fn paf_open(file: &str) -> Result<String> {
    let mut paf_raw: String = String::with_capacity(10_000_000_000);
    let mut paf = File::open(&Path::new(file))?;
    paf.read_to_string(&mut paf_raw)?;
    Ok(paf_raw)
}

// Return (read id of query, start position of query, end position of query, query length) and the
// same information of read target.
pub fn parse<'a>(
    line: &'a str,
) -> Option<(&'a str, usize, usize, usize, &'a str, usize, usize, usize)> {
    let contents: Vec<&str> = line.split('\t').collect();
    let query_length: usize = contents[1].parse().ok()?;
    let query_start: usize = contents[2].parse().ok()?;
    let query_end: usize = contents[3].parse().ok()?;
    let target_length: usize = contents[6].parse().ok()?;
    let target_start: usize = contents[7].parse().ok()?;
    let target_end: usize = contents[8].parse().ok()?;
    let _alignment_block_length: usize = contents[10].parse().ok()?;
    Some((
        &contents[0],
        query_start,
        query_end,
        query_length,
        &contents[5],
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

use std::collections::HashSet;
fn get_ids(file: &str) -> HashSet<String> {
    BufReader::new(std::fs::File::open(&std::path::Path::new(file)).unwrap())
        .lines()
        .skip(1)
        .filter_map(|e| e.ok())
        .collect()
}

pub fn paf_file_to_intervals_with_id(paf: &str, ids: &str) -> Vec<(String, Interval)> {
    let ids = get_ids(ids);
    let mut summary: Summary = HashMap::new();
    let reader = BufReader::new(File::open(&Path::new(paf)).unwrap());
    for line in reader.lines().filter_map(|e| e.ok()) {
        if let Some((read1_id, read1_s, read1_e, length1, read2_id, read2_s, read2_e, length2)) =
            parse(&line)
        {
            if !ids.contains(read1_id) || !ids.contains(read2_id) {
                continue;
            }
            insert_or_update(&mut summary, read1_id, (read1_s, read1_e), length1);
            insert_or_update(&mut summary, read2_id, (read2_s, read2_e), length2);
        }
    }
    eprintln!("Finish read file. {} reads.", summary.len());
    summary
        .into_par_iter()
        .map(|(id, mappings)| (id, Interval::new(&mappings.0, mappings.1)))
        .collect()
}

pub fn paf_file_to_intervals(paf: &str) -> Vec<(String, Interval)> {
    let reader = BufReader::new(File::open(&Path::new(paf)).unwrap());
    to_intervals(reader)
}

type Summary = HashMap<String, (Vec<(usize, usize)>, usize)>;
#[inline]
fn insert_or_update(summary: &mut Summary, read_id: &str, tuple: (usize, usize), len: usize) {
    if summary.contains_key(read_id) {
        summary.get_mut(read_id).unwrap().0.push(tuple);
    } else {
        summary.insert(read_id.to_string(), (vec![tuple], len));
    }
}

fn to_intervals<R: BufRead>(reader: R) -> Vec<(String, Interval)> {
    let mut summary: Summary = HashMap::new();
    for line in reader.lines().filter_map(|e| e.ok()) {
        if let Some((read1_id, read1_s, read1_e, length1, read2_id, read2_s, read2_e, length2)) =
            parse(&line)
        {
            insert_or_update(&mut summary, read1_id, (read1_s, read1_e), length1);
            insert_or_update(&mut summary, read2_id, (read2_s, read2_e), length2);
        }
    }
    summary
        .into_par_iter()
        .map(|(id, mappings)| (id, Interval::new(&mappings.0, mappings.1)))
        .collect()
}

pub fn string_to_intervals(input: &str) -> Vec<(String, Interval)> {
    let reader = BufReader::new(input.as_bytes());
    to_intervals(reader)
}

pub fn stdin_to_intervals() -> Vec<(String, Interval)> {
    let reader = BufReader::new(std::io::stdin());
    to_intervals(reader)
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
