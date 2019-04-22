extern crate rayon;
extern crate select_mitochondrial_genome;
use rayon::prelude::*;
use select_mitochondrial_genome::Interval;
use std::io::Result;
// trim lower and 3 %
//const THRESHOLD: Option<usize> = Some(3);
// Full.
const THRESHOLD: Option<usize> = None;
#[derive(Debug)]
struct Summary {
    mean: f64,
    sd: f64,
    length: usize,
}

impl Summary {
    fn new(sum: usize, sumsq: usize, length: usize) -> Self {
        let mean = sum as f64 / length as f64;
        let variance = sumsq as f64 / length as f64 - mean * mean;
        Summary {
            mean: mean,
            sd: variance.sqrt(),
            length: length,
        }
    }
    fn summing_up(intervals: &Interval) -> (usize, usize) {
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
    fn update(current: usize, change: i8) -> usize {
        if change < 0 {
            current - 1
        } else {
            current + 1
        }
    }
    // Input:intervals Output: sum, sum of sqare, number of considered position.
    fn summing_up_trim(interval: &Interval, is_trim_on: usize) -> (usize, usize, usize) {
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
    fn from_interval(interval: Interval) -> Self {
        if let Some(is_trim_on) = THRESHOLD {
            let (sum, sumsq, length) = Self::summing_up_trim(&interval, is_trim_on);
            Summary::new(sum, sumsq, length)
        } else {
            let (sum, sumsq) = Self::summing_up(&interval);
            Summary::new(sum, sumsq, interval.length())
        }
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let summary_of_each_read = summarize_coverage(&args[1], &args[2]);
    eprintln!("Finish proc. Output.");
    for (id, summary) in summary_of_each_read {
        println!(
            "{}\t{}\t{}\t{}",
            id, summary.mean, summary.sd, summary.length
        );
    }
    eprintln!("Finish output. Exit.");
    Ok(())
}

fn summarize_coverage(paf_file: &str, id_file: &str) -> Vec<(String, Summary)> {
    select_mitochondrial_genome::paf_file_to_intervals_with_id(paf_file, id_file)
        .into_par_iter()
        .map(|(id, interval)| (id, Summary::from_interval(interval)))
        .collect()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    #[test]
    fn test_summing_up() {
        let res = Interval::from_raw(&vec![(0, 1), (10, -1)], 10);
        let (sum, sumsq) = Summary::summing_up(&res);
        assert_eq!((sum, sumsq), (10, 10));
        let res = Interval::from_raw(&vec![(0, 1), (2, -1), (3, 1), (4, -1)], 4);
        let (sum, sumsq) = Summary::summing_up(&res);
        assert_eq!((sum, sumsq), (3, 3));
        let interval = Interval::new(&vec![(0, 3), (1, 10), (2, 5)], 10);
        let (sum, sumsq) = Summary::summing_up(&interval);
        assert_eq!((sum, sumsq), (15, 27));
        let interval = Interval::new(&vec![(0, 1), (2, 3), (4, 5)], 5);
        let (sum, sumsq) = Summary::summing_up(&interval);
        assert_eq!((sum, sumsq), (3, 3));
        let interval = Interval::new(&vec![(0, 10), (1, 8), (3, 7), (5, 6)], 8);
        let (sum, sumsq) = Summary::summing_up(&interval);
        assert_eq!((sum, sumsq), (22, 58));
    }
}
