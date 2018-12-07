extern crate bio;
extern crate rayon;
use bio::io::fasta;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::Path;
const THR:usize = 50;
fn main() -> std::io::Result<()> {
    let args:Vec<_> = std::env::args().collect();
    let file = fasta::Reader::from_file(&Path::new(&args[1]))?;
    let is_rev = (args.len() > 2) && ("rev" == args[2]);
    let record = file.records().filter_map(|e| e.ok()).nth(0).unwrap();
    let tot = record.seq().len();
    let now = std::time::Instant::now();
    let res = fasta_to_consective_runs(record,is_rev);
    let elapsed = std::time::Instant::now();
    let mut stdout = BufWriter::new(std::io::stdout());
    let mut result = String::new();
    for (index, cons) in res {
        for (s, len) in cons {
            if is_rev{
                result.push_str(&format!("{}\t{}\t{},rev\n",index,tot-len-s,len));
            }else{
                result.push_str(&format!("{}\t{}\t{},temp\n",index,s+index,len));
            }
        }
    }
    writeln!(stdout, "{}",result);
    let fin = std::time::Instant::now();
    eprintln!(
        "{:?} = {:?} + {:?}",
        fin - now,
        elapsed - now,
        fin - elapsed
    );
    Ok(())
}


///(starting point of consecutive true, its length)
type Run = (usize, usize);
type Runs = Vec<Run>;

fn fasta_to_consective_runs(record: fasta::Record, is_rev:bool) -> Vec<(usize, Vec<(usize, usize)>)> {
    let seq = record.seq();
    (1..seq.len())
        .into_par_iter()
        .map(|i| (i, consecutive_true(&diagonal_vec(i, &seq, is_rev))))
        .collect::<Vec<_>>()
}

fn diagonal_vec(i: usize, seq: &[u8], is_rev:bool) -> Vec<bool> {
    if !is_rev{
        (0..seq.len() - i).map(|l| seq[l+i] == seq[l]).collect()
    }else{
        let revcomp = bio::alphabets::dna::revcomp(seq);
        (0..seq.len() - i).map(|l| seq[l+i] == revcomp[l]).collect()
    }
}

fn consecutive_true(rest: &[bool]) -> Runs {
    let mut x = None;
    let mut sofar = Vec::with_capacity(rest.len() / 2);
    for (index, &head) in rest.iter().enumerate() {
        x = match x {
            None if head => Some((index, 1)),
            None if !head => None,
            Some((start, len)) if head => Some((start, len + 1)),
            Some((start, len)) if !head => {
                if len > THR{
                    sofar.push((start, len));
                }
                None
            }
            _ => unreachable!(),
        };
    }
    match x {
        Some(res) if res.1 > THR => sofar.push(res),
        _ => {},
    };
    sofar
}

#[test]
fn test() {
    let testvec = [true, true, false, true, true, true];
    assert_eq!(consecutive_true(&testvec), vec![(0, 2), (3, 3)]);
    let testvec = [false, false, false];
    assert_eq!(consecutive_true(&testvec), vec![]);
    let testvec = [true];
    assert_eq!(consecutive_true(&testvec), vec![(0, 1)]);
    let testvec = [true, true];
    assert_eq!(consecutive_true(&testvec), vec![(0, 2)]);
}
