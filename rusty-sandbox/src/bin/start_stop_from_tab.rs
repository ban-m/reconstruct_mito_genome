extern crate bio_utils;
use bio_utils::lasttab::LastTAB;
use std::collections::HashMap;
fn get_first_alignment(tabs: &[LastTAB]) -> Option<(&str, usize)> {
    let tab = tabs.iter().min_by_key(|t| t.seq2_start_from_forward())?;
    let id = tab.seq1_name();
    let position = tab.seq1_start_from_forward();
    Some((id, position))
}

fn get_last_alignment(tabs: &[LastTAB]) -> Option<(&str, usize)> {
    let tab = tabs.iter().min_by_key(|t| t.seq2_end_from_forward())?;
    let id = tab.seq1_name();
    let position = tab.seq1_end_from_forward();
    Some((id, position))
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let tabs: Vec<_> = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter(|line| !line.starts_with("#"))
        .filter_map(|line| LastTAB::from_line(&line))
        .filter(|tab| tab.seq1_matchlen() > 1000 && tab.seq2_matchlen() > 1000)
        .collect();
    eprintln!("{}", tabs.len());
    let references: Vec<_> = {
        let mut reference: Vec<_> = tabs
            .iter()
            .map(|r| {
                let id = r.seq1_name();
                let len = r.seq1_len();
                (id, len)
            })
            .collect();
        reference.sort_by_key(|e| e.1);
        reference.dedup();
        reference
    };
    eprintln!("{:?}", references);
    let mut start_stop_count: HashMap<String, Vec<usize>> = references
        .into_iter()
        .map(|(id, len)| (id.to_string(), vec![0; len + 1]))
        .collect();
    let mut reads: HashMap<_, Vec<_>> = HashMap::new();
    for tab in tabs {
        let entry = reads
            .entry(tab.seq2_name().to_string())
            .or_insert(Vec::new());
        entry.push(tab);
    }
    for (_id, tabs) in reads.into_iter() {
        if let Some((ref_name, position)) = get_first_alignment(&tabs) {
            if let Some(res) = start_stop_count.get_mut(ref_name) {
                res[position] += 1;
            }
        }
        if let Some((ref_name, position)) = get_last_alignment(&tabs) {
            if let Some(res) = start_stop_count.get_mut(ref_name) {
                res[position] += 1;
            }
        }
    }
    for (id, counts) in start_stop_count {
        for (pos, count) in counts.iter().enumerate() {
            println!("{}\t{}\t{}", id, pos, count);
        }
    }
    Ok(())
}
