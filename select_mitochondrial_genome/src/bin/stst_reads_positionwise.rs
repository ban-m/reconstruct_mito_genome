extern crate rayon;
extern crate select_mitochondrial_genome;
use rayon::prelude::*;
use select_mitochondrial_genome::*;
use std::fs::File;
use std::io::Result;
use std::io::{BufWriter, Write};
use std::path::Path;
#[derive(Debug)]
struct StStSummary {
    stop_reads: Vec<(usize, u8)>,
    start_reads: Vec<(usize, u8)>,
}

impl StStSummary {
    fn from_vec_by(interval: &Vec<(usize, i8)>, which: i8, _len: usize) -> Vec<(usize, u8)> {
        let (mut current_pos, mut count, mut counts) = (0, 0, vec![]);
        for &(pos, _) in interval.iter().filter(|(_, ch)| ch == &which) {
            if current_pos == pos {
                count += 1;
            } else{
                if count != 0{
                    counts.push((current_pos, count));
                }
                current_pos = pos;
                count = 1;
            }
        }
        if count != 0{
            counts.push((current_pos, count));
        }
        counts
    }

    fn from_interval_inner(
        interval: Vec<(usize, i8)>,
        len: usize,
    ) -> (Vec<(usize, u8)>, Vec<(usize, u8)>) {
        let starts = Self::from_vec_by(&interval, 1, len);
        let ends = Self::from_vec_by(&interval, -1, len);
        (starts, ends)
    }

    fn from_interval(interval: Interval) -> Self {
        let length = interval.length();
        let interval = interval
            .inner()
            .iter()
            .map(|&(pos, ch)| if ch < 0 { (pos - 1, ch) } else { (pos, ch) })
            .collect();
        let (starts, stops) = Self::from_interval_inner(interval, length);
        StStSummary {
            stop_reads: stops,
            start_reads: starts,
        }
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let coverage_of_each_read = positionwise_stst(&args[1]);
    let mut start = BufWriter::new(File::create(&Path::new(&args[2]))?);
    let mut stop = BufWriter::new(File::create(&Path::new(&args[3]))?);
    for (id, stst) in coverage_of_each_read {
        write!(&mut start, "{}", id)?;
        print_vec(&mut start, &stst.start_reads)?;
        write!(&mut stop, "{}", id)?;
        print_vec(&mut stop, &stst.stop_reads)?;
    }
    Ok(())
}

fn print_vec<W>(wtr: &mut BufWriter<W>, rec: &Vec<(usize, u8)>) -> Result<()>
where
    W: Write,
{
    write!(wtr, ",{}", rec.len())?;
    for (pos, _) in rec {
        write!(wtr, ",{}", pos)?;
    }
    for (_, count) in rec {
        write!(wtr, ",{}", count)?;
    }
    writeln!(wtr, "")
}

fn positionwise_stst(paf: &str) -> Vec<(String, StStSummary)> {
    paf_file_to_intervals(paf)
        .into_par_iter()
        .map(|(id, interval)| (id, StStSummary::from_interval(interval)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn from_vec_by_test() {
        let test = vec![(1, 1), (3, 1), (4, 1), (4, 1), (5, 1)];
        assert_eq!(
            StStSummary::from_vec_by(&test, 1, 6),
            vec![(1,1),(3,1),(4,2),(5,1)]
        );
        assert_eq!(StStSummary::from_vec_by(&test, -1, 6), vec![]);
        let test = vec![(1, -1), (1, -1), (2, -1), (2, -1), (2, -1), (5, -1)];
        assert_eq!(
            StStSummary::from_vec_by(&test, -1, 6),
            vec![(1,2),(2,3),(5,1)]
        );
    }

}
