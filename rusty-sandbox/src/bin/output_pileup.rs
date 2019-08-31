extern crate bio;
extern crate bio_utils;
extern crate rusty_sandbox;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
const NAME_LENGTH: usize = 17; // "[00000 00000Clip]"
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let sam_records = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    let fastq_records = rusty_sandbox::open_fastq_into_hashmap(&args[2])?;
    eprintln!("Opened dataset");
    let output = |target: String| -> std::io::Result<()> {
        let pileup = get_pileup_of(&sam_records[&target], &fastq_records, &target);
        eprintln!("Finish pileup. Pileup is {} height.", pileup.len());
        let mut file = loop {
            print!("Enter file to write:");
            std::io::stdout().flush()?;
            let file = get_input()?;
            match File::create(&Path::new(&file)).map(BufWriter::new) {
                Ok(res) => break res,
                Err(why) => eprintln!("Error:{:?}", why),
            };
        };
        let window = 200;
        let len = pileup[0].len();
        let mut start = 0;
        while start < len {
            let end = (start + window).min(len);
            for line in &pileup {
                writeln!(&mut file, "{}", String::from_utf8_lossy(&line[start..end]))?;
            }
            writeln!(&mut file)?;
            start += window;
        }
        file.flush()
    };
    loop {
        print!("Enter read ID(or \"q\" to exit.):");
        std::io::stdout().flush()?;
        let target = get_input()?;
        if sam_records.contains_key(&target) {
            if let Err(why) = output(target) {
                eprintln!("{:?}", why);
                break;
            }
        } else if target == "q" {
            println!("Quitting...");
            break;
        } else if target == "LONGEST" {
            println!("Searching the longest read...");
            let target = fastq_records
                .iter()
                .filter(|(name, _)| sam_records.contains_key(name.as_str()))
                .max_by_key(|(_, seq)| seq.seq().len())
                .unwrap()
                .0
                .clone();
            println!("target:{}, {}", target, sam_records.contains_key(&target));
            if let Err(why) = output(target) {
                eprintln!("{:?}", why);
                break;
            }
        } else {
            println!("Invalid Read ID");
            println!("The IDs below are possible candidates");
            for key in sam_records.keys().take(10) {
                println!("{}", key);
            }
        }
    }
    Ok(())
}

fn get_input() -> std::io::Result<String> {
    let mut input = String::new();
    std::io::stdin().read_line(&mut input)?;
    Ok(input.trim().to_string())
}

use bio::io::fastq::Record;
use bio_utils::sam;
use std::collections::BinaryHeap;
fn get_pileup_of(
    sams: &[sam::Sam],
    records: &HashMap<String, Record>,
    target: &str,
) -> Vec<Vec<u8>> {
    let coordinate: Vec<_> = format!("{:width$}", target, width = NAME_LENGTH)
        .bytes()
        .chain(records[target].seq().iter().flat_map(|x| vec![*x, b'-']))
        .chain(format!("{:width$}", target, width = NAME_LENGTH).bytes())
        .collect();
    let len = coordinate.len();
    // eprintln!("Coordinates was constructed. {} Length.", len);
    let alignments = concat_split_reads(sams, records);
    // eprintln!("Alignments were recovered.");
    let mut pileups: BinaryHeap<Run> = BinaryHeap::new();
    for rip in alignments {
        let pos = rip.start;
        match pileups.peek() {
            Some(lane) if lane.vacant < pos => {}
            _ => pileups.push(Run::new(pos, len)),
        };
        let mut lane = pileups.pop().unwrap();
        // eprintln!("Overwite {}-{}", pos, pos + rip.alignment.len());
        for (idx, &x) in rip.alignment.iter().enumerate() {
            lane[idx + pos] = x;
        }
        lane.vacant = rip.end;
        pileups.push(lane);
    }
    // eprintln!("Pushed all the alignments.");
    let mut result = pileups.into_vec();
    result.sort_by_key(|e| e.first);
    result.reverse();
    let mut result: Vec<_> = result.into_iter().map(|e| e.content).collect();
    result.push(coordinate);
    result
}

// concat reads in splited to combine.
fn concat_split_reads(sams: &[sam::Sam], records: &HashMap<String, Record>) -> Vec<ReadInPileup> {
    let mut raw: Vec<&sam::Sam> = sams.iter().collect();
    raw.sort_by(|a, b| match (a.q_name()).cmp(b.q_name()) {
        Ordering::Equal => match (a.is_template()).cmp(&b.is_template()) {
            Ordering::Equal => (a.pos()).cmp(&b.pos()),
            x => x,
        },
        x => x,
    });
    // eprintln!("Sorted {} alignments.", raw.len());
    let mut alignments = vec![];
    let mut prev = &sams[0];
    let mut buf: Vec<&sam::Sam> = vec![&sams[0]];
    for sam in &raw[1..] {
        if sam.q_name() == prev.q_name() && sam.is_template() == prev.is_template() {
            buf.push(sam);
        } else {
            alignments.extend(ReadInPileup::from(&buf, records[prev.q_name()].seq()));
            buf.clear();
            buf.push(sam);
            prev = sam;
        }
    }
    alignments.extend(ReadInPileup::from(&buf, records[prev.q_name()].seq()));
    alignments.sort_by_key(|e| e.start);
    alignments
}

#[derive(Debug, Clone)]
struct ReadInPileup {
    // Including header and footer.
    alignment: Vec<u8>,
    start: usize,
    end: usize,
}

impl ReadInPileup {
    fn from(sams: &[&sam::Sam], seq: &[u8]) -> Vec<Self> {
        let name = sams[0].q_name();
        // eprintln!("Summing up {}'s alignment({} in total)", name, sams.len());
        let seq = if sams[0].is_template() {
            seq.to_vec()
        } else {
            bio::alphabets::dna::revcomp(seq)
        };
        let mut mappings: BinaryHeap<Map> = BinaryHeap::new();
        for mapping in sams.iter().map(|sam| Map::from_sam(sam, &seq)) {
            let mergable = match mappings.peek() {
                Some(e) => e.mergable(&mapping),
                _ => false,
            };
            if mergable {
                let mut m = mappings.pop().unwrap();
                m.merge(mapping);
                mappings.push(m);
            } else {
                mappings.push(mapping);
            }
        }
        let mut res = vec![];
        while let Some(mapping) = mappings.pop() {
            res.push(ReadInPileup::from_map(mapping, name));
        }
        res
    }
    fn from_map(map: Map, name: &str) -> Self {
        let start = 2 * map.start;
        let end = start + map.alignment.len() + NAME_LENGTH;
        let alignment: Vec<_> = format!("[{:05} {:05}Clip]", name, map.head_clip)
            .bytes()
            .chain(map.alignment)
            .chain(format!("[{:05} {:05}Clip]", name, map.tail_clip).bytes())
            .collect();
        ReadInPileup {
            alignment,
            start,
            end,
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
struct Map {
    head_clip: usize,
    tail_clip: usize,
    // This is the position at the reference, not the coordinate in the output.
    // In short, it is just sam.pos().
    start: usize,
    end: usize,
    alignment: Vec<u8>,
}

impl Map {
    fn from_sam(sam: &sam::Sam, seq: &[u8]) -> Self {
        use sam::Op;
        let cigar: Vec<sam::Op> = sam.cigar();
        let head_clip = match cigar.first() {
            Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
            _ => 0,
        };
        let tail_clip = match cigar.last() {
            Some(Op::SoftClip(l)) | Some(Op::HardClip(l)) => *l,
            _ => 0,
        };
        let (start, end) = sam.get_range();
        let mut alignment: Vec<u8> = vec![];
        // cursour
        let mut idx = head_clip;
        for op in cigar {
            match op {
                Op::Align(l) | Op::Match(l) | Op::Mismatch(l) => {
                    for x in &seq[idx..idx + l] {
                        alignment.push(*x);
                        alignment.push(b' ');
                    }
                    idx += l;
                }
                Op::Insertion(l) => {
                    alignment.pop();
                    alignment.push(b'0' + l.min(9) as u8);
                    idx += l;
                }
                Op::Deletion(l) => {
                    for _ in &seq[idx..idx + l] {
                        alignment.push(b'-');
                        alignment.push(b' ');
                    }
                }
                _ => {}
            }
        }
        Map {
            head_clip,
            tail_clip,
            start,
            end,
            alignment,
        }
    }
    fn mergable(&self, other: &Self) -> bool {
        self.end < other.start && self.head_clip + self.alignment.len() / 2 > other.head_clip
    }
    fn merge(&mut self, other: Self) {
        self.tail_clip = other.tail_clip;
        self.end = other.end;
        self.alignment.push(b' ');
        for _ in self.end..(other.start - 2) {
            self.alignment.push(b'-');
            self.alignment.push(b' ');
        }
    }
}

impl Ord for Map {
    fn cmp(&self, other: &Self) -> Ordering {
        other.end.cmp(&self.end)
    }
}

impl PartialOrd for Map {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct Run {
    content: Vec<u8>,
    // the location at which the last character inserted.
    vacant: usize,
    // the location at which the first character inserted.
    first: usize,
}

impl Ord for Run {
    fn cmp(&self, other: &Self) -> Ordering {
        other.vacant.cmp(&self.vacant)
    }
}
impl PartialOrd for Run {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Run {
    fn new(first: usize, len: usize) -> Self {
        Run {
            content: vec![b' '; len],
            first,
            vacant: first,
        }
    }
}

impl std::ops::Index<usize> for Run {
    type Output = u8;
    fn index(&self, index: usize) -> &Self::Output {
        &self.content[index]
    }
}

impl std::ops::IndexMut<usize> for Run {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.content[index]
    }
}
