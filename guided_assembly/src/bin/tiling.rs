// Tiling
// Create tiling from LAST's TAB format output
// The output format is a line-delimited format, with each line consisting of
// [read id]:([assign])*( [assign])+, where [read id] is a unique identifier of each read(i.e., fastq id),
// [assign] is one of [contig name]-[bin index] or G-[gap length].
// For example,
// m1213/23243/43423/0_2324:0-1 0-2 0-3 G-2343 2-0 2-1 2-2

//(score, ctg name, ctg start, ctg aln len, read id, read start, read aln len)
struct LastAln<'a>{
    score:usize,
    ctgname:usize,
    ctgstart:usize,
    ctgend:usize,
    read_id:&'a str,
    read_start:usize,
    read_end:usize,
}

const BINWIDTH: usize = 200;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::{BufWriter, Write};
use std::path::Path;

enum Assignment {
    Tile(usize, usize),
    Gap(usize),
}

impl std::fmt::Display for Assignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Assignment::Tile(ctg, num) => write!(f, "{}-{}", ctg, num),
            Assignment::Gap(size) => write!(f, "G-{}", size),
        }
    }
}

fn open_file(file: &str) -> std::io::Result<String> {
    let mut input = String::new();
    let mut file = File::open(&Path::new(file))?;
    file.read_to_string(&mut input)?;
    Ok(input)
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut wtr = BufWriter::new(std::io::stdout());
    for (read_id, assignments) in tiling(open_file(&args[1])?) {
        write!(&mut wtr, "{}:", read_id)?;
        if assignments.is_empty(){
            writeln!(&mut wtr,"");
        }else{
            for i in 0..assignments.len() - 1 {
                write!(&mut wtr, "{} ", assignments[i])?;
            }
            writeln!(&mut wtr, "{}", assignments.last().unwrap())?;
        }
    }
    Ok(())
}

fn aggregate<'a> (input: &'a str) -> HashMap<String,Vec<LastAln<'a>>>{
    let mut alignments:HashMap<String,Vec<_>> = HashMap::new();
    for (id, aln) in input.lines().filter(|e| !e.starts_with('#'))
        .map(|e|e.split('\t').collect())
        .filter_map(|e|parse(e))
        .map(|e|(e.read_id.to_string(),e)){
            eprint!("{}",&id);
            let entry = alignments.entry(id).or_default();
            if entry.is_empty(){
                eprintln!("-> First in");
            }else{
                // eprintln!("Already has {} alignments",entry.len());
            }
            entry.push(aln);
        }
    eprintln!("Parsing Finish");
    alignments
}

fn tiling(input: String) -> Vec<(String, Vec<Assignment>)> {
    aggregate(&input).into_iter()
        .filter_map(|(read_id,aligns)|Some((read_id,trim_overlap(aligns)?)))
        .map(|(read_id,aligns)|(read_id, tile_read(aligns)))
        .collect()
}

// (score, contig name, start, alnlen, strand, length, readid, readstart, readalnlen, strand, readlen)
fn parse<'a>(contents: Vec<&'a str>) -> Option<LastAln<'a>> {
    if contents.len() < 10 {
        eprintln!("{:?} Please check", &contents);
        return None
    }
    let parse_return = |x:&str| x.parse::<usize>().ok();
    let score:usize = parse_return(contents[0])?;
    let ctgname: usize = parse_return(contents[1])?;
    let ctgstart: usize = parse_return(contents[2])?;
    let ctg_aln_length: usize = parse_return(contents[3])?;
    if ctg_aln_length < BINWIDTH {
        return None;
    }
    // To tile the read, one need to scale the alignmnet size;
    let _ctg_strand = if contents[4] == "+" { '+' } else { '-'};
    let _ctglen = parse_return(contents[5])?;
    let read_id = &contents[6];
    let read_start = parse_return(contents[7])?;
    let read_aln_length  = parse_return(contents[8])?;
    let read_strand = if contents[9] == "+" { '+' }else{ '-'};
    let read_length = parse_return(contents[10])?;
    if read_strand == '+' {
        Some(LastAln{
            score:score,
            ctgname:ctgname,
            ctgstart:ctgstart,
            ctgend:ctgstart + ctg_aln_length,
            read_id:read_id,
            read_start:read_start,
            read_end:read_start + read_aln_length,
        })
    }else{
        Some(LastAln{
            score:score,
            ctgstart:ctgstart,
            ctgname:ctgname,
            ctgend:ctgstart + ctg_aln_length,
            read_id:read_id,
            read_start:read_length - read_start -read_aln_length,
            read_end:read_length - read_start
        })
    }
}
fn trim_overlap(mut alignment:Vec<LastAln>) -> Option<Vec<LastAln>> {
    alignment.sort_by_key(|e|e.read_start);
    let mut no_overlap = vec![];
    let mut alignment = alignment.into_iter().peekable();
    let mut current  = alignment.next()?;
    while let Some(next_aln) = alignment.peek(){
        if next_aln.read_start < current.read_end {
            // Overlap detected. Choose one.
            if next_aln.score < current.score {
                // retain current one. consume next one and go on to next iteration.
                alignment.next();
            }else{
                // Discard current alignmnet and go on to next iteration.
                current = alignment.next().unwrap();
            }
        }else{
            // No overlap. Retain.
            no_overlap.push(current);
            current = alignment.next().unwrap();
        }
    }
    Some(no_overlap)
}

fn tile_read(aligns:Vec<LastAln>)->Vec<Assignment>{
    let mut result = vec![];
    aligns.windows(2)
        .for_each(|e|{
            let prev = &e[0];
            let next = &e[1];
            if next.read_start - prev.read_end > BINWIDTH{
                result.push(Assignment::Gap(next.read_start - prev.read_end));
            }
            result.extend(into_assignment(prev));
        });
    result
}

fn into_assignment(align:&LastAln) -> Vec<Assignment>{
    (align.ctgstart / BINWIDTH .. align.ctgend / BINWIDTH)
        .map(|bin|
             Assignment::Tile(align.ctgname,bin))
        .collect()
}



