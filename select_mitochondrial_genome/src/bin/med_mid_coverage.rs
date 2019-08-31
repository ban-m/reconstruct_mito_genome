extern crate rayon;
extern crate select_mitochondrial_genome;
use rayon::prelude::*;
use select_mitochondrial_genome::Interval;
use std::collections::HashMap;
use std::io::Result;

#[derive(Debug)]
struct Summary {
    median: usize,
    mode: usize,
    length: usize,
}

impl Summary {
    fn new(median: usize, mode: usize, length: usize) -> Self {
        Summary {
            median,
            mode,
            length,
        }
    }
    #[inline]
    fn update(current: usize, change: i8) -> usize {
        if change < 0 {
            current - 1
        } else {
            current + 1
        }
    }
    // Input:interval, Return (median,mode)
    fn summing_up(interval: &Interval) -> (usize, usize) {
        let (mut current_position, mut coverage) = (0, 0);
        let mut coverages = vec![];
        for &(pos, coverage_change) in interval.inner().iter() {
            if pos > current_position {
                (current_position..pos).for_each(|_| coverages.push(coverage));
                current_position = pos;
            }
            coverage = Self::update(coverage, coverage_change);
        }
        (current_position..interval.length()).for_each(|_| coverages.push(coverage));
        coverages.sort();
        (median(&coverages, interval.length()), mode(&coverages))
    }
    fn from_interval(interval: &Interval) -> Self {
        let (median, mode) = Self::summing_up(interval);
        Summary::new(median, mode, interval.length())
    }
}

// Input should be sorted.
fn median(cov: &Vec<usize>, length: usize) -> usize {
    cov[length / 2]
}

// Input should be sorted.
fn mode(cov: &Vec<usize>) -> usize {
    let (mut current, mut num) = (0, 0);
    let (mut most_occured, mut max_num) = (0, 0);
    for x in cov {
        if &current == x {
            num += 1;
        } else {
            if max_num < num {
                max_num = num;
                most_occured = current;
            }
            current = *x;
            num = 0;
        }
    }
    if max_num < num {
        max_num = num;
        most_occured = current;
    }
    debug_assert!(max_num != 0);
    most_occured
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let summary_of_each_read = summarize_coverage(&args[1]);
    for (id, summary) in summary_of_each_read {
        println!(
            "{}\t{}\t{}\t{}",
            id, summary.median, summary.mode, summary.length
        );
    }
    Ok(())
}

fn summarize_coverage(paf: &str) -> HashMap<String, Summary> {
    select_mitochondrial_genome::paf_file_to_intervals(paf)
        .into_par_iter()
        .map(|(id, interval)| (id, Summary::from_interval(&interval)))
        .collect()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    #[test]
    fn mode_test() {
        let input = vec![1, 1, 1, 1, 1];
        assert_eq!(mode(&input), 1);
        let input = vec![1, 2, 2, 3];
        assert_eq!(mode(&input), 2);
        let input = vec![1, 2, 2, 3, 3, 3, 4];
        assert_eq!(mode(&input), 3);
        let input = vec![1, 2, 2, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10];
        assert_eq!(mode(&input), 10);
    }
    #[test]
    fn median_test() {
        let input = vec![0, 1, 2, 3, 4, 5];
        assert_eq!(median(&input, 6), 3);
        let input = vec![21, 21, 24];
        assert_eq!(median(&input, 3), 21);
        let input: Vec<usize> = (0..100).collect();
        assert_eq!(median(&input, 100), 50);
    }
}
