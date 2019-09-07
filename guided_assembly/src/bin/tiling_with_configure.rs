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
    ctgname: String,
    ctgstart: usize,
    ctgend: usize,
    read_id: String,
    read_start: usize,
    read_end: usize,
    is_forward: bool,
}

const BINWIDTH: usize = 200;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

enum Assignment {
    Tile(usize, u8, usize),
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

fn parse_last_aln(contents: String) -> Option<LastAln> {
    if contents.split('\t').count() < 10 {
        panic!("{:?} Please check", &contents);
    }
    let parse_return = |x: &str| x.parse::<usize>().ok();
    let contents: Vec<_> = contents.split('\t').collect();
    let score: usize = parse_return(contents[0])?;
    let ctgname: String = contents[1].to_string();
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
            is_forward: true,
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
            is_forward: false,
        })
    }
}

fn parse_peak_file(file: &str) -> std::io::Result<(HashMap<String, usize>, Vec<Vec<(usize, u8)>>)> {
    let mut contents: Vec<(String, usize)> = File::open(file)
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let mut e = e.split('\t');
            Some((e.next()?.to_owned(), e.next().and_then(|e| e.parse().ok())?))
        })
        .collect();
    contents.sort();
    let index = {
        let mut num = 0;
        contents.iter().fold(HashMap::new(), |mut res, (name, _)| {
            if !res.contains_key(name) {
                res.insert(name.to_string(), num);
                num += 1;
            }
            res
        })
    };
    eprintln!("{:?}", &index);
    let mut peaks: Vec<Vec<(usize, u8)>> = vec![vec![]; index.len()];
    for (name, position) in contents {
        let idx = index[&name];
        let color = if peaks[idx].is_empty() {
            0
        } else {
            peaks[idx].last().unwrap().1 + 1
        };
        eprintln!(
            "Fill {}'s {}-{} by {}",
            name,
            peaks[idx].len(),
            position,
            color
        );
        peaks[idx].push((position, color));
    }
    for peak in peaks.iter_mut() {
        peak.sort_by_key(|e| e.0);
    }
    eprintln!("{:?}", peaks);
    Ok((index, peaks))
}

fn trim_overlap(mut alns: Vec<LastAln>) -> Option<Vec<LastAln>> {
    alns.sort_by_key(|e| e.read_start);
    let mut no_overlap = vec![];
    let mut alns = alns.into_iter().peekable();
    let mut current = alns.next()?;
    while let Some(next) = alns.peek() {
        if next.read_start < current.read_end {
            // Overlap detected. Choose one.
            if next.score < current.score {
                // retain current one. consume next one and go on to next iteration.
                alns.next();
            } else {
                // Discard current alignmnet and go on to next iteration.
                current = alns.next().unwrap();
            }
        } else {
            // No overlap. Retain.
            no_overlap.push(current);
            current = alns.next().unwrap();
        }
    }
    Some(no_overlap)
}

fn tile_read(
    alns: Vec<LastAln>,
    index: &HashMap<String, usize>,
    peak: &[Vec<(usize, u8)>],
) -> Vec<Assignment> {
    let mut result = vec![];
    alns.windows(2).for_each(|e| {
        let prev = &e[0];
        let next = &e[1];
        if next.read_start - prev.read_end > BINWIDTH {
            result.push(Assignment::Gap(next.read_start - prev.read_end));
        }
        let contig_idx = index[&prev.ctgname];
        let peak = &peak[contig_idx];
        result.extend(into_assignment(prev, contig_idx, peak));
    });
    result
}

fn into_assignment(align: &LastAln, index: usize, peak: &[(usize, u8)]) -> Vec<Assignment> {
    if align.is_forward {
        (align.ctgstart / BINWIDTH..align.ctgend / BINWIDTH)
            .map(|bin| {
                let position = bin * BINWIDTH;
                let (p, c) = match peak.binary_search_by_key(&position, |e| e.0) {
                    Ok(res) => peak[res],
                    Err(res) => peak[res - 1],
                };
                let bin = (position - p) / BINWIDTH;
                Assignment::Tile(index, c, bin)
            })
            .collect()
    } else {
        (align.ctgstart / BINWIDTH..align.ctgend / BINWIDTH)
            .rev()
            .map(|bin| {
                let position = bin * BINWIDTH;
                let (p, c) = match peak.binary_search_by_key(&position, |e| e.0) {
                    Ok(res) => peak[res],
                    Err(res) => peak[res - 1],
                };
                let bin = (position - p) / BINWIDTH;
                Assignment::Tile(index, c, bin)
            })
            .collect()
    }
}

fn aggregate(alignments: Vec<LastAln>) -> HashMap<String, Vec<LastAln>> {
    let mut result: HashMap<String, Vec<_>> = HashMap::new();
    for aln in alignments {
        let entory = result.entry(aln.read_id.clone()).or_default();
        entory.push(aln);
    }
    result
}

fn tiling(
    index: &HashMap<String, usize>,
    alignments: Vec<LastAln>,
    peaks: &Vec<Vec<(usize, u8)>>,
) -> Vec<(String, Vec<Assignment>)> {
    let chunks = aggregate(alignments);
    chunks
        .into_iter()
        .filter_map(|(read_id, aligns)| Some((read_id, trim_overlap(aligns)?)))
        .map(|(read_id, aligns)| (read_id, tile_read(aligns, &index, &peaks)))
        .collect()
}

// First argument is the mappings(LAST format), and the second is peak call format.
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (index, peaks) = parse_peak_file(&args[2])?;
    eprintln!("Parsed");
    eprintln!("{:?},{}", &index, peaks.len());
    let alignments = File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
        .filter(|e| !e.starts_with('#'))
        .filter_map(parse_last_aln)
        .collect();
    let out = std::io::stdout();
    let mut out = BufWriter::new(out.lock());
    let tiles = tiling(&index, alignments, &peaks);
    eprintln!("Tiled");
    for (read_id, assignments) in tiles {
        let mut record = read_id + ":";
        for assignment in assignments {
            record.push_str(&format!("{} ", assignment));
        }
        record.push('\n');
        out.write_all(record.as_bytes())?;
    }
    Ok(())
}
