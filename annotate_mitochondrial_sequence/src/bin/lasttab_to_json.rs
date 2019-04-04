extern crate annotate_mitochondrial_sequence;
extern crate bio;
extern crate bio_types;
extern crate rayon;
extern crate serde;
extern crate serde_json;
use annotate_mitochondrial_sequence::{open_annotation, Annotation, Map};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let last = trim_overlap(open_last(&args[1])?).unwrap();
    let annotation: HashMap<_, _> = open_annotation(&args[2])?;
    let maps: Vec<_> = last
        .into_par_iter()
        .filter_map(|aln| convert_to_map(&aln, &annotation))
        .collect();
    println!("{}", serde_json::ser::to_string_pretty(&maps).unwrap());
    Ok(())
}

fn overlap(start1: usize, end1: usize, start2: usize, end2: usize) -> usize {
    if end1 <= start2 || end2 <= start1 {
        // [s1,e1) and [s2,e2) does not overlap.
        return 0;
    }
    assert!(end1.min(end2) >= start1.max(start2));
    end1.min(end2) - start1.max(start2)
}

fn calc_confidence(aln: &LastAln, annot: &Annotation) -> f64 {
    let gene_length = (annot.end - annot.start) as f64;
    overlap(annot.start, annot.end, aln.ref_start, aln.ref_end) as f64 / gene_length
}

fn find_largest_ovlp_from_bucket<'a>(
    aln: &LastAln,
    bucket: &'a Vec<Annotation>,
) -> Option<(&'a Annotation, f64)> {
    let index = match bucket.binary_search_by_key(&aln.ref_start, |e| e.start) {
        Ok(res) => res,
        Err(res) if res < bucket.len() => res,
        Err(_) => return None,
    };
    let (mut max_confidence, mut largest_overlap) = (-100.0, &bucket[index]);
    for annot in bucket
        .iter()
        .skip(index)
        .take_while(|annot| annot.start < aln.ref_end)
    {
        let confidence = calc_confidence(aln, annot);
        if confidence > max_confidence {
            max_confidence = confidence;
            largest_overlap = annot;
        }
    }
    if max_confidence == -100.0 {
        // can not find any overlaps.
        None
    } else {
        Some((largest_overlap, max_confidence))
    }
}

fn find_largest_overlap<'a>(
    aln: &LastAln,
    annot: &'a HashMap<String, Vec<Annotation>>,
) -> Option<(&'a Annotation, f64)> {
    let bucket = annot.get(&aln.ref_name)?;
    find_largest_ovlp_from_bucket(aln, &bucket)
}

fn convert_to_map(aln: &LastAln, annot: &HashMap<String, Vec<Annotation>>) -> Option<Map> {
    let (largest_ovlp, confidence) = find_largest_overlap(aln, annot)?;
    let is_forward = largest_ovlp.strand * aln.strand == 1;
    let (start,end) = (aln.query_start.min(aln.query_end),
                       aln.query_end.max(aln.query_start));
    Some(Map::new(
        &aln.query_name,
        false,
        is_forward,
        start as u64,
        end as u64,
        &aln.ref_name,
        &largest_ovlp.gene,
        confidence,
    ))
}

struct LastAln {
    score: usize,
    ref_name: String,
    ref_start: usize,
    ref_end: usize,
    query_name: String,
    query_start: usize,
    query_end: usize,
    strand: i8,
}

fn open_last(file: &str) -> std::io::Result<Vec<LastAln>> {
    let input = open_file(file)?;
    let res: Vec<_> = input
        .lines()
        .filter(|e| !e.starts_with('#'))
        .filter_map(|e| {
            let e: Vec<_> = e.split('\t').collect();
            LastAln::parse(e)
        })
        .collect();
    Ok(res)
}

fn open_file(file: &str) -> std::io::Result<String> {
    let mut input = String::new();
    let mut file = File::open(&Path::new(file))?;
    file.read_to_string(&mut input)?;
    Ok(input)
}

impl LastAln {
    fn parse(contents: Vec<&str>) -> Option<Self> {
        if contents.len() < 10 {
            eprintln!("{:?} Please check", &contents);
            return None;
        }
        let parse_return = |x: &str| x.parse::<usize>().ok();
        let score: usize = parse_return(contents[0])?;
        let ctgname = contents[1].to_string();
        let ctgstart: usize = parse_return(contents[2])?;
        let ctg_aln_length: usize = parse_return(contents[3])?;
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
                ref_name: ctgname,
                ref_start: ctgstart,
                ref_end: ctgstart + ctg_aln_length,
                query_name: read_id,
                query_start: read_start,
                strand:1,
                query_end: read_start + read_aln_length,
            })
        } else {
            Some(LastAln {
                score: score,
                ref_start: ctgstart,
                ref_name: ctgname,
                ref_end: ctgstart + ctg_aln_length,
                query_name: read_id,
                query_start: read_length - read_start - read_aln_length,
                strand:-1,
                query_end: read_start - read_start,
            })
        }
    }
}

fn trim_overlap(mut alignment: Vec<LastAln>) -> Option<Vec<LastAln>> {
    if alignment.is_empty(){
        eprintln!("Error: the LAST alignment file is empty. Please check the input. Return");
        return None;
    }
    alignment.sort_by_key(|e| e.query_start);
    let mut no_overlap = vec![];
    let mut alignment = alignment.into_iter().peekable();
    let mut current = alignment.next()?;
    while let Some(next_aln) = alignment.peek() {
        if next_aln.query_start < current.query_end {
            // Overlap detected. Choose one.
            if next_aln.score < current.score {
                // retain current one. consume next one and go on to next iteration.
                alignment.next();
            } else {
                // Discard current alignmnet and go on to next iteration.
                current = alignment.next().unwrap();
            }
        } else {
            // No overlap. Retain.
            no_overlap.push(current);
            current = alignment.next().unwrap();
        }
    }
    Some(no_overlap)
}
