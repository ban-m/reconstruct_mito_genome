// Tiling
// Create tiling from LAST's TAB format output
// The output format is a line-delimited format, with each line consisting of
// [read id]:([assign])*( [assign])+, where [read id] is a unique identifier of each read(i.e., fastq id),
// [assign] is one of [contig index]-[peak index]-[bin index] or G-[gap length].
// For example,
// m1213/23243/43423/0_2324:0-0-1 0-0-2 0-0-3 G-2343 2-1-0 2-0-1 2-0-2

//(score, ctg name, ctg start, ctg aln len, read id, read start, read aln len)
struct LastAln {
    score: usize,
    ctgname: usize,
    ctgstart: usize,
    ctgend: usize,
    read_id: String,
    read_start: usize,
    read_end: usize,
}

const BINWIDTH: usize = 200;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

enum Assignment {
    Tile(usize, usize, usize),
    Gap(usize),
}

impl std::fmt::Display for Assignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Assignment::Tile(ctg, peak, num) => write!(f, "{}-{}-{}", ctg, peak, num),
            Assignment::Gap(size) => write!(f, "G-{}", size),
        }
    }
}

fn parse_last_aln(contents: String, _index: HashMap<String, usize>) -> Option<LastAln> {
    if contents.len() < 10 {
        eprintln!("{:?} Please check", &contents);
        return None;
    }
    let parse_return = |x: &str| x.parse::<usize>().ok();
    let contents: Vec<_> = contents.split('\t').collect();
    let score: usize = parse_return(contents[0])?;
    let ctgname: usize = parse_return(contents[1])?;
    let ctgstart: usize = parse_return(contents[2])?;
    let ctg_aln_length: usize = parse_return(contents[3])?;
    if ctg_aln_length < BINWIDTH {
        return None;
    }
    // To tile the read, one need to scale the alignmnet size;
    let _ctg_strand = if contents[4] == "+" { '+' } else { '-' };
    let _ctglen = parse_return(contents[5])?;
    let read_id = contents[6].to_string();
    let read_start = parse_return(contents[7])?;
    let read_aln_length = parse_return(contents[8])?;
    let read_strand = if contents[9] == "+" { '+' } else { '-' };
    let read_length = parse_return(contents[10])?;
    if read_strand == '+' {
        Some(LastAln {
            score: score,
            ctgname: ctgname,
            ctgstart: ctgstart,
            ctgend: ctgstart + ctg_aln_length,
            read_id: read_id,
            read_start: read_start,
            read_end: read_start + read_aln_length,
        })
    } else {
        Some(LastAln {
            score: score,
            ctgstart: ctgstart,
            ctgname: ctgname,
            ctgend: ctgstart + ctg_aln_length,
            read_id: read_id,
            read_start: read_length - read_start - read_aln_length,
            read_end: read_length - read_start,
        })
    }
}

fn parse_peak_file(file: &str) -> std::io::Result<(HashMap<String, usize>, Vec<Vec<u8>>)> {
    let contents: Vec<String> = File::open(file)
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .collect();
    let index = contents
        .iter()
        .skip(1)
        .filter_map(|e| e.split('\t').nth(0))
        .enumerate()
        .fold(HashMap::new(), |mut res, (idx, name)| {
            if !res.contains_key(name) {
                res.insert(name.to_string(), idx);
            }
            res
        });
    let mut peaks: Vec<_> = vec![vec![]; index.len()];
    for line in contents {
        let mut e = line.split('\t');
        let name = e.next().unwrap();
        let position: usize = e.next().and_then(|e| e.parse().ok()).unwrap();
        let idx = index[name];
        let color = if peaks[idx].is_empty() {
            0
        } else {
            *peaks[idx].last().unwrap()
        };
        eprintln!(
            "Fill {}'s {}-{} by {}",
            name,
            peaks[idx].len(),
            position,
            color
        );
        (peaks[idx].len()..position).for_each(|e| peaks[idx].push(color));
    }
    Ok((index, peaks))
}




// First argument is the mappings(LAST format), and the second is peak call format.
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (index, peaks) = parse_peak_file(&args[2])?;
    let alignments = File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| parse_last_aln(e, index));
    let (read_id, assignments) in tiling(index, alignments,peaks);
    Ok(())
}
