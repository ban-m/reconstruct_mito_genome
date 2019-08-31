extern crate bio;
extern crate bio_utils;

use bio_utils::sam::Sam;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
const MIN_CLIP: usize = 100;
fn open_sam_with_contig_name(path: &str, contig_name: &str) -> std::io::Result<Vec<Sam>> {
    let mut reader = BufReader::new(std::fs::File::open(&Path::new(path))?).lines();
    let (first_line, contigs) = {
        let mut header = HashSet::new();
        loop {
            let line = reader.next().unwrap().unwrap();
            if line.starts_with("@SQ") {
                let name = line.split('\t').nth(1).unwrap().split(':').nth(1).unwrap();
                header.insert(name.to_string());
            } else if !line.starts_with("@") {
                break (line, header);
            }
        }
    };
    if !contigs.contains(contig_name) {
        eprintln!("Contig name:{} does not exists", contig_name);
        eprintln!("The available contig keys as follows:");
        for key in &contigs {
            eprintln!("{}", key);
        }
        eprintln!("Exit.");
        return Err(std::io::Error::from(std::io::ErrorKind::Other));
    }
    let sam_records: Vec<_> = vec![first_line]
        .into_iter()
        .chain(reader.filter_map(|e| e.ok()))
        .filter_map(|e| bio_utils::sam::Sam::new(&e))
        .filter(|e| e.r_name() == contig_name)
        .collect();
    let primary_queries: HashSet<_> = sam_records
        .iter()
        .filter(|e| e.is_primary())
        .map(|e| e.q_name().to_string())
        .collect();
    let mut sam_records: Vec<_> = sam_records
        .into_iter()
        .filter(|e| primary_queries.contains(e.q_name()))
        .collect();
    sam_records.sort_by(|a, b| a.q_name().cmp(b.q_name()));
    Ok(sam_records)
}

fn coverage(
    sam_records: &Vec<bio_utils::sam::Sam>,
    reference: &bio::io::fasta::Record,
) -> Vec<u16> {
    let mut coverage: Vec<u16> = vec![0; reference.seq().len()];
    for record in sam_records {
        let mut pos = record.pos() - 1;
        use bio_utils::sam::Op;
        for op in record.cigar_ref() {
            match &&op {
                Op::Align(l) | Op::Mismatch(l) | Op::Match(l) => {
                    for i in pos..(pos + l) {
                        coverage[i] += 1;
                    }
                    pos += l;
                }
                Op::Deletion(l) => pos += l,
                _ => {}
            }
        }
    }
    coverage
}

fn extract_clippings(sam: &bio_utils::sam::Sam) -> Vec<usize> {
    use bio_utils::sam::Op;
    let cigar_ref = sam.cigar_ref();
    let head = match cigar_ref.first() {
        Some(Op::HardClip(l)) | Some(Op::SoftClip(l)) => *l,
        _ => 0,
    };
    let tail = match cigar_ref.last() {
        Some(Op::HardClip(l)) | Some(Op::SoftClip(l)) => *l,
        _ => 0,
    };
    vec![head, tail]
}

fn determine_clip_type(
    records: &[&bio_utils::sam::Sam],
    threshold: usize,
    query_length: usize,
) -> (usize, usize, usize) {
    let mut range: Vec<(usize, usize)> = records
        .iter()
        .filter_map(|e| {
            let (start, end) = e.mapped_region();
            // Remove proper alignmnet.
            if start < threshold && query_length - end < threshold {
                None
            } else {
                if e.is_template() {
                    Some((start, end))
                } else {
                    Some((query_length - end, query_length - start))
                }
            }
        })
        .collect();
    range.sort_by_key(|e| e.0);
    let check = |(start, end)| {
        if start > threshold && query_length - end > threshold {
            (0, 1, 0)
        } else {
            (1, 0, 0)
        }
    };
    if range.is_empty() {
        (0, 0, 0)
    } else if range.len() == 1 {
        check(range[0])
    } else {
        let coverage = calc_coverage(&range);
        if coverage > query_length.max(2 * threshold) - 2 * threshold {
            (0, 0, 1)
        } else {
            range
                .into_iter()
                .map(check)
                .fold((0, 0, 0), |(a, b, c), (x, y, z)| (a + x, b + y, c + z))
        }
    }
}

fn calc_coverage(ranges: &[(usize, usize)]) -> usize {
    let mut cumsum = 0;
    let (mut p_s, mut p_e) = ranges[0];
    for (start, end) in &ranges[1..] {
        if end < &p_e {
            // contained.
        } else if start < &p_e {
            // Overlapping elongation.
            p_e = *end
        } else {
            // Pop.
            cumsum += p_e - p_s;
            p_s = *start;
            p_e = *end;
        }
    }
    cumsum + (p_e - p_s)
}

fn clippings(
    sam_records: &Vec<bio_utils::sam::Sam>,
    queries: &HashMap<String, bio::io::fastq::Record>,
    _reference: &bio::io::fasta::Record,
) -> (usize, usize, usize, usize, Vec<usize>) {
    // Note that `sam_records` are already sorted array.
    let clipping_length: Vec<_> = sam_records.iter().flat_map(extract_clippings).collect();
    let threshold = clipping_length[clipping_length.len() * 75 / 100].max(MIN_CLIP);
    let mut previous = "";
    let mut buf: Vec<&bio_utils::sam::Sam> = vec![];
    let mut query_length = 0;
    let (mut head_tail, mut double_clip, mut junction) = (0, 0, 0);
    for record in sam_records.iter() {
        if previous == record.q_name() {
            buf.push(record);
        } else {
            if !buf.is_empty() {
                let (h_t, d_c, junc) = determine_clip_type(&buf, threshold, query_length);
                head_tail += h_t;
                double_clip += d_c;
                junction += junc;
                buf.clear();
            }
            previous = record.q_name();
            buf.push(record);
            query_length = queries[previous].seq().len();
        }
    }
    let (h_t, d_c, junc) = determine_clip_type(&buf, threshold, query_length);
    head_tail += h_t;
    double_clip += d_c;
    junction += junc;
    (head_tail, double_clip, junction, threshold, clipping_length)
}

#[inline]
fn is_different_base(x: u8, y: u8) -> bool {
    match (x, y) {
        (b'a', b'a')
        | (b'A', b'a')
        | (b'a', b'A')
        | (b'A', b'A')
        | (b'c', b'c')
        | (b'C', b'c')
        | (b'c', b'C')
        | (b'C', b'C')
        | (b'g', b'g')
        | (b'G', b'g')
        | (b'g', b'G')
        | (b'G', b'G')
        | (b't', b't')
        | (b'T', b't')
        | (b't', b'T')
        | (b'T', b'T') => false,
        _ => true,
    }
}

fn error_profile(
    sam_records: &Vec<bio_utils::sam::Sam>,
    queries: &HashMap<String, bio::io::fastq::Record>,
    reference: &bio::io::fasta::Record,
) -> Vec<(usize, usize, usize, usize)> {
    let reference = reference.seq();
    let mut error_profile: Vec<(usize, usize, usize, usize)> = vec![(0, 0, 0, 0); reference.len()];
    for record in sam_records {
        let query = if record.is_template() {
            queries[record.q_name()].seq().to_vec()
        } else {
            bio::alphabets::dna::revcomp(queries[record.q_name()].seq())
        };
        let mut ref_pos = record.pos() - 1;
        let mut query_pos = 0;
        use bio_utils::sam::Op;
        for &op in record.cigar_ref() {
            match op {
                Op::Align(l) | Op::Match(l) | Op::Mismatch(l) => {
                    for i in 0..l {
                        if is_different_base(query[query_pos + i], reference[ref_pos + i]) {
                            error_profile[ref_pos + i].0 += 1;
                        }
                        error_profile[ref_pos + i].3 += 1;
                    }
                    query_pos += l;
                    ref_pos += l;
                }
                Op::Deletion(l) => {
                    for i in 0..l {
                        error_profile[ref_pos + i].1 += 1;
                        error_profile[ref_pos + i].3 += 1;
                    }
                    ref_pos += l;
                }
                Op::Insertion(l) => {
                    error_profile[ref_pos].2 += l;
                    query_pos += l;
                }
                Op::SoftClip(l) | Op::HardClip(l) => query_pos += l,
                _ => {}
            };
        }
    }
    error_profile
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        eprintln!("[sam file] [query file] [reference file] [output directly]  [contig name]");
        return Ok(());
    }
    let sam_records = open_sam_with_contig_name(&args[1], &args[5])?;
    println!("Sam records:{} records", sam_records.len());
    let good_id: HashSet<_> = sam_records.iter().map(|e| e.q_name().to_string()).collect();
    let mut out = PathBuf::from(&args[4]);
    if !out.exists() {
        std::fs::create_dir_all(&out)?;
    }
    let queries: HashMap<_, _> = bio::io::fastq::Reader::from_file(&Path::new(&args[2]))?
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| good_id.contains(e.id()))
        .map(|e| (e.id().to_string(), e))
        .collect();
    println!("There are {} queries.", queries.len());
    let reference = bio::io::fasta::Reader::from_file(&Path::new(&args[3]))?
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| e.id() == &args[5])
        .nth(0)
        .unwrap();
    // Coverage calculate
    let coverage = coverage(&sam_records, &reference);
    {
        let zeros = coverage.iter().filter(|&&e| e == 0).count();
        println!(
            "{}/{}={}",
            zeros,
            reference.seq().len(),
            zeros as f64 / reference.seq().len() as f64
        );
        out.push("coverage.tsv");
        let mut wtr = BufWriter::new(std::fs::File::create(&out)?);
        for (idx, x) in coverage.iter().enumerate() {
            writeln!(&mut wtr, "{}\t{}", idx, x)?;
        }
        let sum: u64 = coverage.iter().map(|&e| e as u64).sum();
        let sumsq: u64 = coverage.iter().map(|&e| (e as u64).pow(2)).sum();
        let mean: f64 = sum as f64 / reference.seq().len() as f64;
        let var: f64 = sumsq as f64 / reference.seq().len() as f64 - mean * mean;
        println!("Mean:{}\tVar:{}", mean, var);
        out.pop();
    }
    // Clippings
    let (head_tail, double_clip, junctions, threshold, clip_lengths) =
        clippings(&sam_records, &queries, &reference);
    {
        out.push("clippings.dat");
        println!("Clipping threshold:{}", threshold);
        println!(
            "HeadTail:{}\nDouble:{}\nJunctions:{}",
            head_tail, double_clip, junctions
        );
        let mut wtr = BufWriter::new(std::fs::File::create(&out)?);
        for x in &clip_lengths {
            writeln!(&mut wtr, "{}", x)?;
        }
        out.pop();
    }
    // Error Profiles
    let error_profile = error_profile(&sam_records, &queries, &reference);
    {
        out.push("error_profile.dat");
        let mut wtr = BufWriter::new(std::fs::File::create(&out)?);
        for (idx, (subst, ins, del, depth)) in error_profile.into_iter().enumerate() {
            writeln!(&mut wtr, "{}\t{}\t{}\t{}\t{}", idx, subst, ins, del, depth)?;
        }
    }
    Ok(())
}
