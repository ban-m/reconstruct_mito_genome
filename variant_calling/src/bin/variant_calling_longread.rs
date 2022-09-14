extern crate bio_utils;
extern crate clap;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate env_logger;
use clap::{App, Arg};
use std::fs::File;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let matches = App::new("naive_VC_longreads")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Call variant from pileups")
        .arg(
            Arg::with_name("reads")
                .required(true)
                .short('r')
                .long("reads")
                .value_name("READS")
                .help("Raw long reads<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("reference")
                .required(true)
                .short('l')
                .long("reference")
                .value_name("REFERENCE")
                .help("Reference sequence<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("alignments")
                .required(true)
                .short('b')
                .long("bam")
                .value_name("BAM")
                .help("alignments file<BAM>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("all")
                .short('a')
                .long("all")
                .help("Dump all locations on the reference."),
        )
        .arg(
            Arg::with_name("freq_thr")
                .short('f')
                .long("freq_thr")
                .takes_value(true)
                .default_value(&"0.7")
                .help("Frequency threshold at each position"),
        )
        .arg(
            Arg::with_name("depth")
                .short('d')
                .long("depth")
                .takes_value(true)
                .default_value(&"30")
                .help("Sequence depth threshould at each position."),
        )
        .arg(
            Arg::with_name("gff")
                .short('g')
                .long("gff")
                .takes_value(true)
                .help("GFF file to add genes containing variants in the output.)"),
        )
        .arg(
            Arg::with_name("gff_convert")
                .short('n')
                .long("name")
                .takes_value(true)
                .help("A tab-delimited file containig <Fasta ID>-><GFF ID>"),
        )
        .arg(
            Arg::with_name("verbose")
                .short('v')
                .multiple(true)
                .help("Output debug to the standard error."),
        )
        .get_matches();
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let reads = matches
        .value_of("reads")
        .map(|file| match bio_utils::fasta::parse_into_vec(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let reference = matches
        .value_of("reference")
        .map(|file| match bio_utils::fasta::parse_into_vec(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let alignments = matches
        .value_of("alignments")
        .map(|file| match rust_htslib::bam::Reader::from_path(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let freq_thr: f64 = matches
        .value_of("freq_thr")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let depth: u32 = matches
        .value_of("depth")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let is_all = matches.is_present("all");
    let gff = {
        if let Some(mut gff) = matches.value_of("gff").map(|file| match open_gff(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        }) {
            debug!("GFFfile opened");
            if let Some(file) = matches.value_of("gff_convert") {
                debug!("GFF genome name would be converted by {}", file);
                let mapping = match open_gff_convert(file) {
                    Ok(res) => res,
                    Err(why) => panic!("{}:{}", why, file),
                };
                for (key, value) in mapping.iter() {
                    debug!("{}->{}", key, value);
                }
                for record in gff.iter_mut() {
                    if let Some(new_name) = mapping.get(&record.0) {
                        record.0 = new_name.clone();
                    }
                }
            }
            Some(gff)
        } else {
            None
        }
    };
    debug!("FreqThr:{}", freq_thr);
    debug!("Depth Thre:{}", depth);
    let variants = variant_calling(&reads, &reference, alignments, is_all, (freq_thr, depth));
    use std::io::{BufWriter, Write};
    let wtr = std::io::stdout();
    let mut wtr = BufWriter::new(wtr.lock());
    for (refname, vars) in variants {
        for (pos, count, pileup, refr, query) in vars {
            let (r, q) = (refr as char, query as char);
            let mut line = format!("{}\t{}\t{}\t{}\t{}\t{}", refname, pos, pileup, count, r, q);
            if let Some(ref gff) = gff {
                line.push_str("\t");
                match find(gff, &refname, pos) {
                    Some(desc) => line.push_str(desc),
                    None => line.push_str("none"),
                };
            }
            writeln!(&mut wtr, "{}", line)?;
        }
    }
    Ok(())
}

fn find<'a>(gff: &'a [GFF], refname: &str, pos: usize) -> Option<&'a str> {
    gff.iter()
        .filter_map(|&(ref name, start, end, ref desc)| {
            if name == refname && start <= pos && end <= pos {
                Some(desc.as_str())
            } else {
                None
            }
        })
        .next()
}

use bio_utils::fasta::Record;
use rust_htslib::bam::Read;
use rust_htslib::bam::Reader;
type VariantCallResult = Vec<(String, Vec<(usize, u32, u32, u8, u8)>)>;
fn variant_calling(
    reads: &[Record],
    reference: &[Record],
    mut bam: Reader,
    is_all: bool,
    param: (f64, u32),
) -> VariantCallResult {
    let header = bam.header().clone();
    let reads: HashMap<_, _> = reads
        .iter()
        .map(|e| (e.id().as_bytes().to_vec(), e))
        .collect();
    // A,C,G,T, Ins, Del,
    let mut counts: Vec<_> = (0..header.target_count())
        .filter_map(|tid| header.target_len(tid))
        .map(|len| vec![[0; 6]; len as usize + 1])
        .collect();
    for r in bam
        .records()
        .filter_map(|e| e.ok())
        .filter(|r| !r.is_secondary() && r.mapq() > 20)
    {
        let refr = &mut counts[r.tid() as usize];
        let read = if r.is_reverse() {
            revcmp(reads[r.qname()].seq())
        } else {
            reads[r.qname()].seq().to_vec()
        };
        let mut rpos = r.pos() as usize;
        let mut qpos = 0 as usize;
        use rust_htslib::bam::record::Cigar;
        for cigar in r.cigar().iter() {
            match cigar {
                Cigar::Del(d) => {
                    for i in (0..*d).map(|i| i as usize) {
                        refr[rpos + i][5] += 1;
                    }
                    rpos += *d as usize;
                }
                Cigar::Diff(m) | Cigar::Equal(m) | Cigar::Match(m) => {
                    for i in (0..*m).map(|i| i as usize) {
                        match read[qpos + i] {
                            b'A' | b'a' => refr[rpos + i][0] += 1,
                            b'C' | b'c' => refr[rpos + i][1] += 1,
                            b'G' | b'g' => refr[rpos + i][2] += 1,
                            b'T' | b't' => refr[rpos + i][3] += 1,
                            _ => unreachable!(),
                        }
                    }
                    qpos += *m as usize;
                    rpos += *m as usize;
                }
                Cigar::HardClip(c) | Cigar::SoftClip(c) => {
                    qpos += *c as usize;
                }
                Cigar::Ins(c) => {
                    for i in (0..*c).map(|i| i as usize) {
                        refr[rpos + i][4] += 1;
                    }
                    qpos += *c as usize;
                }
                Cigar::Pad(_) | Cigar::RefSkip(_) => {}
            }
        }
    }
    counts
        .into_iter()
        .enumerate()
        .filter(|x| !(x.1).is_empty())
        .map(|(tid, vs)| (tid as u32, vs))
        .filter_map(|(tid, vs)| Some((String::from_utf8(header.tid2name(tid).to_vec()).ok()?, vs)))
        .filter_map(|(refname, counts)| {
            let reference = reference.iter().find(|r| refname == r.id())?;
            // assert_eq!(reference.seq().len(), counts.len());
            let vs = variant_calling_filter(reference.seq(), &counts, param, is_all);
            Some((refname, vs))
        })
        .collect()
}

fn variant_calling_filter(
    reference: &[u8],
    counts: &[[u32; 6]],
    (freq_thr, depth_thr): (f64, u32),
    is_all: bool,
) -> Vec<(usize, u32, u32, u8, u8)> {
    let mut called = vec![];
    for (pos, (&refbase, count)) in reference.iter().zip(counts).enumerate() {
        let (argmax, max) = get_max_base(&count[..4]);
        let is_diff_base = refbase.to_ascii_uppercase() != argmax;
        let depth = count.iter().take(4).sum::<u32>();
        let is_freq_large = max as f64 / depth as f64 > freq_thr;
        let is_depth_enough = depth > depth_thr;
        if (is_diff_base && is_freq_large && is_depth_enough) || is_all {
            called.push((pos, depth, max, refbase, argmax));
        }
        let insertion = count[4];
        let is_insertion = insertion as f64 / depth as f64 > freq_thr;
        if is_insertion && is_depth_enough {
            called.push((pos, depth, insertion, b'-', b'I'));
        }
        let deletion = count[5];
        let is_deletion = deletion as f64 / depth as f64 > freq_thr;
        if is_deletion && is_depth_enough {
            called.push((pos, depth, deletion, b'-', b'D'));
        }
    }
    called
}
fn get_max_base(bc: &[u32]) -> (u8, u32) {
    let (argmax, &max) = bc.iter().enumerate().max_by_key(|e| e.1).unwrap();
    match argmax {
        0 => (b'A', max),
        1 => (b'C', max),
        2 => (b'G', max),
        3 => (b'T', max),
        _ => unreachable!(),
    }
}

#[inline]
fn revcmp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&e| match e {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

type GFF = (String, usize, usize, String);
fn open_gff(file: &str) -> std::io::Result<Vec<GFF>> {
    let res: Vec<_> = File::open(file)
        .map(BufReader::new)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|line| {
            let line: Vec<_> = line.split('\t').collect();
            if line.len() < 5 {
                None
            } else {
                let chr = line[0].to_string();
                let start: usize = line[3].parse().ok()?;
                let end: usize = line[4].parse().ok()?;
                let desc = line.last()?.to_string();
                Some((chr, start, end, desc))
            }
        })
        .collect();
    Ok(res)
}
use std::collections::HashMap;
fn open_gff_convert(file: &str) -> std::io::Result<HashMap<String, String>> {
    File::open(file).map(|f| {
        BufReader::new(f)
            .lines()
            .filter_map(|e| e.ok())
            .map(|line| {
                let line: Vec<_> = line.split('\t').collect();
                (line[0].to_string(), line[1].to_string())
            })
            .collect()
    })
}
