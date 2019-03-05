extern crate rayon;
extern crate rust_htslib;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::io::Result;
use std::path::Path;
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let coverages: Vec<(_, _)> = positionwise_coverage(&args[1])?;
    for (id, coverage) in coverages {
        let csv = vec_to_string(coverage);
        println!("{},{}", id, csv);
    }
    Ok(())
}

#[derive(Debug)]
struct Alignment {
    seq1: u32,
    seq2: u32,
    seq1_start: usize,
    seq1_end: usize,
    seq2_start: usize,
    seq2_end: usize,
}

fn positionwise_coverage(file: &str) -> Result<Vec<(String, Vec<u32>)>> {
    let mut record = bam::Record::new();
    let mut reader = open_bam(file)?;
    let (name_to_tid, tid_to_name, tid_to_length) = get_maps(&reader);
    let mut matrix: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    while let Ok(_) = reader.read(&mut record) {
        if !record.is_unmapped() {
            let align = to_alignment(&record, &name_to_tid);
            let to_seq1 = matrix.entry(align.seq1).or_insert(vec![]);
            to_seq1.push((align.seq1_start, align.seq1_end));
            let to_seq2 = matrix.entry(align.seq2).or_insert(vec![]);
            to_seq2.push((align.seq2_start, align.seq2_end));
        }
    }
    eprintln!("{:?}", matrix.len());
    let result: Vec<(String, Vec<u32>)> = matrix
        .into_iter()
        .map(|(tid, maps)| {
            eprintln!("{}", tid);
            get_coverage_of(maps, tid_to_length[&tid], &tid_to_name[&tid])
        })
        .collect();
    Ok(result)
}
fn get_maps(
    reader: &bam::Reader,
) -> (
    HashMap<Vec<u8>, u32>,
    HashMap<u32, String>,
    HashMap<u32, usize>,
) {
    let header = reader.header();
    let mut name_to_tid = HashMap::new();
    let mut tid_to_name = HashMap::new();
    let mut tid_to_length = HashMap::new();
    for name in header.target_names() {
        let tid = header.tid(name).unwrap();
        let length = header.target_len(tid).unwrap();
        name_to_tid.insert(name.to_vec(), tid);
        let name = String::from_utf8_lossy(name).to_string();
        eprintln!("{}\t{}\t{}", name, tid, length);
        tid_to_name.insert(tid, name);
        tid_to_length.insert(tid, length as usize);
    }
    (name_to_tid, tid_to_name, tid_to_length)
}

fn to_alignment(record: &bam::Record, name_to_tid: &HashMap<Vec<u8>, u32>) -> Alignment {
    let seq1 = name_to_tid[record.qname()]; // query
    let seq2 = record.tid() as u32; // target
    let (seq1_start, seq1_end, seq2_start, seq2_end) = sum_up_cigar(record);
    // eprintln!(
//        "{}\t{}-{}, {}\t{}-{}",
    //seq1, seq1_start, seq1_end, seq2, seq2_start, seq2_end
    //    );
    Alignment {
        seq1,
        seq2,
        seq1_start,
        seq1_end,
        seq2_start,
        seq2_end,
    }
}

fn sum_up_cigar(record: &bam::Record) -> (usize, usize, usize, usize) {
    let pos = record.pos() as u32;
    use rust_htslib::bam::record::Cigar;
    let (mut seq1_start, mut seq1_length, mut seq2_length) = (0, 0, 0);
    let mut first_clip = true;
    for e in record.cigar().iter() {
        match e {
            Cigar::SoftClip(l) | Cigar::HardClip(l) if first_clip => {
                seq1_start += l;
                first_clip = false;
            },
            Cigar::Del(l) => seq2_length += l,
            Cigar::Ins(l) => seq1_length += l,
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                seq1_length += l;
                seq2_length += l;
            }
            _ => {}
        }
    }
    (
        seq1_start as usize,
        (seq1_start + seq1_length) as usize,
        pos as usize,
        (pos + seq2_length) as usize,
    )
}

fn open_bam(file: &str) -> Result<bam::Reader> {
    match bam::Reader::from_path(&Path::new(file)) {
        Ok(res) => Ok(res),
        Err(why) => {
            eprintln!("{:?}", why);
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "open failed",
            ))
        }
    }
}

fn get_coverage_of(maps: Vec<(usize, usize)>, len: usize, name: &str) -> (String, Vec<u32>) {
    let mut coverage = vec![0; len];
    for (start, end) in maps {
        eprintln!("{},{}-{}", len, start, end);
        for i in start..end {
            coverage[i] += 1;
        }
    }
    (name.to_string(), coverage)
}

fn vec_to_string(coverage: Vec<u32>) -> String {
    let mut res = String::new();
    for (index, val) in coverage.into_iter().enumerate() {
        if index != 0 {
            res.push(',');
        }
        res.push_str(&format!("{}", val));
    }
    res
}
