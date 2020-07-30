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
        .filter_map(|aln| convert_to_maps(&aln, &annotation))
        .flat_map(|aln| aln)
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

fn find_ovlpping_genes_from_species<'a>(
    aln: &'a LastAln,
    genes: &'a [Annotation],
) -> Vec<&'a Annotation> {
    let index = match genes.binary_search_by_key(&aln.ref_start, |e| e.start) {
        Ok(res) => res,
        Err(res) => res,
    };
    genes
        .iter()
        .skip(index)
        .take_while(|annot| annot.start < aln.ref_end)
        .collect()
}

fn find_overlapping_genes<'a>(
    aln: &'a LastAln,
    annot: &'a HashMap<String, Vec<Annotation>>,
) -> Option<Vec<&'a Annotation>> {
    let bucket = annot.get(&aln.ref_name)?;
    Some(find_ovlpping_genes_from_species(aln, &bucket))
}

fn locate_gene_to_contig(aln: &LastAln, annot: &Annotation) -> (u64, u64) {
    let start = aln.get_corresponding_point_of_query(annot.start) as u64;
    let end = aln.get_corresponding_point_of_query(annot.end) as u64;
    eprintln!(
        "Alignment [{},{}) -> [{},{}). Gene:[{},{}).",
        aln.ref_start, aln.ref_end, aln.query_start, aln.query_end, annot.start, annot.end
    );
    eprintln!(
        "Mapping [{},{}) ({} % confidence)",
        start,
        end,
        calc_confidence(aln, annot) * 100.
    );
    (start, end)
}

fn convert_to_map(aln: &LastAln, annot: &Annotation) -> Map {
    let (start, end) = locate_gene_to_contig(aln, annot);
    let confidence = calc_confidence(aln, annot);
    Map::new(
        &aln.query_name,
        false,
        annot.strand * aln.strand == 1,
        start,
        end,
        &aln.ref_name,
        &annot.gene,
        confidence,
    )
}

// Remark: it returns an array of Map, wrapped by option.
// It is because a mapping can be broken into maps of genes.
fn convert_to_maps(aln: &LastAln, annot: &HashMap<String, Vec<Annotation>>) -> Option<Vec<Map>> {
    let result = find_overlapping_genes(aln, annot)?
        .into_iter()
        .map(|ovlp| convert_to_map(aln, ovlp))
        .collect();
    Some(result) // could be empty vector.
}

#[derive(Debug)]
enum Edit {
    Match(usize),
    // (ref,query). Either is zero.
    Diff(usize, usize),
}

impl Edit {
    fn new(e: &str) -> Option<Self> {
        if e.contains(':') {
            let e: Vec<_> = e
                .split(':')
                .filter_map(|e| e.parse::<usize>().ok())
                .collect();
            if e.len() == 2 {
                return Some(Edit::Diff(e[0], e[1]));
            }
        }
        Some(Edit::Match(e.parse::<usize>().ok()?))
    }
}

#[derive(Debug)]
struct LastAln {
    score: usize,
    ref_name: String,
    ref_start: usize,
    ref_end: usize,
    query_name: String,
    query_start: usize,
    query_end: usize,
    strand: i8,
    blocks: Vec<Edit>,
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
    fn parse_blocks(blocks: &str) -> Option<Vec<Edit>> {
        let mut res = vec![];
        for e in blocks.trim().split(',') {
            res.push(Edit::new(e)?);
        }
        Some(res)
    }
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
        let blocks = Self::parse_blocks(contents[11])?;
        if read_strand == '+' {
            Some(LastAln {
                score,
                ref_name: ctgname,
                ref_start: ctgstart,
                ref_end: ctgstart + ctg_aln_length,
                query_name: read_id,
                query_start: read_start,
                query_end: read_start + read_aln_length,
                strand: 1,
                blocks,
            })
        } else {
            Some(LastAln {
                score,
                ref_start: ctgstart,
                ref_name: ctgname,
                ref_end: ctgstart + ctg_aln_length,
                query_name: read_id,
                query_start: read_length - read_start - read_aln_length,
                query_end: read_length - read_start,
                strand: -1,
                blocks: blocks.into_iter().rev().collect(),
            })
        }
    }
    fn get_corresponding_point_of_query(&self, ref_point: usize) -> usize {
        let (mut query_pos, mut ref_pos) = (self.query_start, self.ref_start);
        for ed in self.blocks.iter() {
            match ed {
                Edit::Match(l) => {
                    if ref_pos + l >= ref_point {
                        return query_pos + ref_point - ref_pos;
                    } else {
                        ref_pos += l;
                        query_pos += l;
                    }
                }
                Edit::Diff(r, q) if q == &0 => {
                    if ref_pos + r >= ref_point {
                        return query_pos;
                    } else {
                        ref_pos += r;
                    }
                }
                Edit::Diff(_r, q) => query_pos += q,
            }
        }
        query_pos
    }
}

fn trim_overlap(mut alignment: Vec<LastAln>) -> Option<Vec<LastAln>> {
    if alignment.is_empty() {
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
    for aln in &no_overlap {
        eprintln!("[{},{})", aln.query_start, aln.query_end);
    }
    Some(no_overlap)
}
