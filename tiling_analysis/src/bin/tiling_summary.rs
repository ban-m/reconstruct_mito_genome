extern crate rayon;
use rayon::prelude::*;

fn convert_inner(encodes: &str) -> (Vec<usize>, Vec<(usize, usize)>) {
    let (mut gaps, mut units) = (vec![], vec![]);
    for unit in encodes.split_whitespace() {
        let mut unit = unit.split('-');
        let unit_or_gap = unit.next().unwrap();
        let size = unit.next().and_then(|e| e.parse().ok()).unwrap();
        match unit_or_gap.parse() {
            Err(_) => gaps.push(size),
            Ok(res) => units.push((res, size)),
        }
    }
    (gaps, units)
}

fn convert(line: &str) -> Option<(Vec<usize>, Vec<(usize, usize)>)> {
    let encodes: &str = line.split(':').nth(1)?;
    Some(convert_inner(encodes))
}

fn count_units(units: Vec<(usize, usize)>) -> Vec<(usize, usize, usize)> {
    if units.is_empty() {
        return vec![];
    }
    let (mut current, mut count) = (&units[0], 0);
    let mut res = vec![];
    for unit in &units[1..] {
        if unit == current {
            count += 1;
        } else {
            res.push((current.0, current.1, count));
            current = unit;
            count = 0;
        }
    }
    res.push((current.0, current.1, count));
    res
}

fn summarize_reads(input: String) -> (Vec<usize>, Vec<(usize, usize, usize)>) {
    let lines: Vec<&str> = input.lines().collect();
    let (gaps, mut units) = lines.into_par_iter().filter_map(convert).reduce(
        || (vec![], vec![]),
        |(mut gaps, mut units), (gap, unit)| {
            gaps.extend(gap);
            units.extend(unit);
            (gaps, units)
        },
    );
    units.sort();
    (gaps, count_units(units))
}

fn summarize(file: &str) -> std::io::Result<(Vec<usize>, Vec<(usize, usize, usize)>)> {
    use std::io::Read;
    let mut file = std::fs::File::open(&std::path::Path::new(file))?;
    let mut input = String::new();
    file.read_to_string(&mut input)?;
    Ok(summarize_reads(input))
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (gaps, units) = summarize(&args[1])?;
    println!("Type\tUnit\tCount");
    for gapsize in gaps {
        println!("GAP\t0\t{}", gapsize); // 0 is a dummy field.
    }
    for (unit, subnum, count) in units {
        println!("{}\t{}\t{}", unit, subnum, count);
    }
    Ok(())
}
