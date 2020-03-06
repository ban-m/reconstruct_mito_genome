type Sequence = (String, [usize; 3]);
type PAF = (Sequence, char, Sequence, u32, u32, u32);
use std::collections::HashMap;

fn parse_sequence(chunks: &[&str]) -> Sequence {
    let id = chunks[0].to_string();
    let lens: Vec<usize> = chunks[1..].iter().filter_map(|e| e.parse().ok()).collect();
    let lens = [lens[0], lens[1], lens[2]];
    (id, lens)
}

fn parse(line: String) -> Option<PAF> {
    let line: Vec<_> = line.split('\t').collect();
    let query = parse_sequence(&line[0..4]);
    let strand = line[4].chars().nth(0)?;
    let reference = parse_sequence(&line[5..9]);
    let match_len: u32 = line[9].parse().ok()?;
    let aln_len: u32 = line[10].parse().ok()?;
    let map_qual: u32 = line[11].parse().ok()?;
    Some((query, strand, reference, match_len, aln_len, map_qual))
}

fn get_first_alignment(paf: &[PAF]) -> Option<(&str, usize)> {
    let paf = paf.iter().min_by_key(|p| ((p.0).1)[1])?;
    let id = &(paf.2).0;
    let position = ((paf.2).1)[1];
    Some((id, position))
}

fn get_last_alignment(paf: &[PAF]) -> Option<(&str, usize)> {
    let paf = paf.iter().min_by_key(|p| ((p.0).1)[2])?;
    let id = &(paf.2).0;
    let position = ((paf.2).1)[2];
    Some((id, position))
}
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let paf: Vec<_> = BufReader::new(File::open(&args[1])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(parse)
        .collect();
    let references: Vec<_> = {
        let mut reference: Vec<_> = paf
            .iter()
            .map(|r| {
                let id = &(r.2).0;
                let len = &((r.2).1)[0];
                (id, len)
            })
            .collect();
        reference.sort_by_key(|e| e.1);
        reference.dedup();
        reference
    };
    eprintln!("{:?}", references);
    let mut start_stop_count: HashMap<_, _> = references
        .into_iter()
        .map(|(id, &len)| (id.clone(), vec![0; len + 1]))
        .collect();
    let mut reads: HashMap<_, Vec<_>> = HashMap::new();
    for record in paf {
        let entry = reads.entry(((record.0).0).clone()).or_insert(Vec::new());
        entry.push(record);
    }
    for (_id, pafs) in reads.into_iter() {
        if let Some((ref_name, position)) = get_first_alignment(&pafs) {
            if let Some(res) = start_stop_count.get_mut(ref_name) {
                res[position] += 1;
            }
        }
        if let Some((ref_name, position)) = get_last_alignment(&pafs) {
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
