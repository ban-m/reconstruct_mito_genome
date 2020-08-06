extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
extern crate serde_json;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate mito_assembler;
use bio_utils::fasta;
use clap::{App, Arg, SubCommand};
use mito_assembler::dump_viewer;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
fn subcommand_create_viewer() -> App<'static, 'static> {
    SubCommand::with_name("create_viewer")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Create Viewer for long reads.")
        .arg(
            Arg::with_name("reads")
                .required(true)
                .short("r")
                .long("reads")
                .value_name("READS")
                .help("Raw long reads<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("contigs")
                .required(true)
                .short("c")
                .long("contigs")
                .value_name("CONTIGS")
                .help("Assembled contigs<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("reference")
                .required(true)
                .short("f")
                .long("reference")
                .value_name("REFERENCE")
                .help("The original reference<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("read_alignments")
                .required(true)
                .short("a")
                .long("read_aln")
                .value_name("ALIGNMENT(Read->C)")
                .help("Alignment from reads to contigs<LastTAB>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("contig_alignments")
                .required(true)
                .short("l")
                .long("contig_aln")
                .value_name("ALIGNMENT(C->Ref)")
                .help("Alignments from contigs to the reference<LastTAB>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("output_dir")
                .required(true)
                .short("o")
                .long("output_dir")
                .value_name("OUT DIR")
                .help("Output directry.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("assignments")
                .required(true)
                .short("t")
                .long("assignments")
                .value_name("ASSIGNMENTS")
                .help("Assignmnet of each read<TSV>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_align_length")
                .long("min_align_length")
                .default_value(&"5000")
                .value_name("MINIMUM LENGTH")
                .help("Minimum alignment length")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Output debug to the standard error."),
        )
}

fn subcommand_decompose() -> App<'static, 'static> {
    SubCommand::with_name("decompose")
        .version("0.1")
        .about("To Decompose long reads")
        .arg(
            Arg::with_name("reads")
                .required(true)
                .short("r")
                .long("reads")
                .value_name("READS")
                .help("Raw long reads<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("reference")
                .required(true)
                .short("c")
                .long("contigs")
                .value_name("CONTIGS")
                .help("Reference<FASTA>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("alignments")
                .required(true)
                .short("a")
                .long("alignments")
                .value_name("ALIGNMENTS")
                .help("Alignments between reads and the reference<LAST TAB>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("selfalignments")
                .required(true)
                .short("s")
                .long("self_alignments")
                .value_name("SELF-ALIGNMENTS")
                .help("Self-vs-Self alignemnts of the reference<LAST TAB>")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("outdir")
                .short("o")
                .long("output")
                .required(true)
                .value_name("OUTPUT_DIRECTORY")
                .help("Output directory")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cluster_num")
                .short("c")
                .long("cluster_num")
                .required(false)
                .value_name("CLUSTER_NUM")
                .help("Minimum cluster number.")
                .default_value(&"2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Output debug to the standard error."),
        )
        .arg(
            Arg::with_name("limit")
                .short("l")
                .long("limit")
                .required(false)
                .value_name("LIMIT")
                .help("Maximum Execution time(sec)")
                .default_value(&"72000")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("no_merge")
                .long("no_merge")
                .required(false)
                .help("Do not exec decompose. Just detects SVs."),
        )
        .arg(
            Arg::with_name("resume")
                .long("resume")
                .required(false)
                .value_name("LOG FILE")
                .help("Log files <JSON>")
                .takes_value(true),
        )
}

fn decompose(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("MMMM started. Debug mode.");
    trace!("Tracing is on.");
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
        .map(|file| match last_tiling::parse_tab_file(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let self_aln = matches
        .value_of("selfalignments")
        .map(|file| match last_tiling::parse_tab_file(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let output_dir = matches
        .value_of("outdir")
        .expect("please specify output directry.");
    debug!("All files opened.");
    let config = last_decompose::error_profile::summarize_tab(&alignments, &reads, &reference);
    let cluster_num: usize = matches
        .value_of("cluster_num")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let limit: u64 = matches
        .value_of("limit")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let no_merge = matches.is_present("no_merge");
    let contigs = last_tiling::contig::Contigs::new(reference);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let initial_clusters = last_decompose::initial_clusters(&encoded_reads, &contigs);
    debug!("Initial clusters constructed");
    {
        if let Err(why) = std::fs::create_dir_all(output_dir) {
            error!("Error Occured while outputing reads.");
            error!("{:?}", why);
            error!("This program did not work successfully.");
            error!("Shutting down...");
            std::process::exit(1);
        }
        let start_stop = last_tiling::get_start_stop(&encoded_reads, &contigs);
        let file = format!("{}/start_stop.tsv", output_dir);
        let mut writer = BufWriter::new(std::fs::File::create(&file)?);
        for (contig, stst) in start_stop {
            for (pos, count) in stst {
                writeln!(&mut writer, "{}\t{}\t{}", contig, pos, count)?;
            }
        }
        let file = format!("{}/split_count.tsv", output_dir);
        let mut writer = BufWriter::new(std::fs::File::create(&file)?);
        for cluster in initial_clusters.iter() {
            for cp in cluster.members.iter().filter_map(|e| e.cr.contig_pair()) {
                let count = cp.reads().len();
                let (s1, e1) = cp.contig1().range();
                let (s2, e2) = cp.contig2().range();
                writeln!(&mut writer, "{}\t{}\t{}\t{}\t{}", s1, e1, s2, e2, count)?;
            }
        }
        let dir = format!("{}/viewer", output_dir);
        if let Err(why) = std::fs::create_dir_all(&dir) {
            error!("Error Occured while outputing reads.");
            error!("{:?}", why);
            error!("This program did not work successfully.");
            error!("Shutting down...");
            std::process::exit(1);
        }
        let repeats = last_tiling::into_repeats(&self_aln, &contigs);
        let file = format!("{}/repeats.json", dir);
        let mut writer = BufWriter::new(std::fs::File::create(&file)?);
        let repeats = serde_json::ser::to_string(&repeats).unwrap();
        writeln!(&mut writer, "{}", repeats)?;
        let file = format!("{}/circos.html", dir);
        let mut writer = BufWriter::new(std::fs::File::create(&file)?);
        writeln!(&mut writer, "{}", mito_assembler::template::TEMPLATE)?;
        let file = format!("{}/style.css", dir);
        let mut writer = BufWriter::new(std::fs::File::create(&file)?);
        writeln!(&mut writer, "{}", mito_assembler::template::STYLE)?;
        debug!("Dumped files");
    }
    for c in &initial_clusters {
        let counts = c.ids().len();
        debug!("{:?} having {} reads.", c, counts);
    }
    let cl = cluster_num;
    debug!("Profiled Error Rates:{}", config);
    let results: HashMap<String, u8> = if !no_merge {
        let result = if matches.is_present("resume") {
            let chunked_reads: Vec<last_decompose::assemble::ChunkedRead> = matches
                .value_of("resume")
                .and_then(|e| std::fs::File::open(e).map(BufReader::new).ok())
                .and_then(|e| serde_json::de::from_reader(e).ok())
                .unwrap();
            last_decompose::assemble::assemble_reads(chunked_reads)
        } else {
            let result = last_decompose::decompose(
                encoded_reads,
                &initial_clusters,
                &contigs,
                &config,
                cl,
                limit,
            );
            let filename = format!("{}/encoded_reads.json", output_dir);
            if let Ok(mut wtr) = std::fs::File::create(&filename).map(std::io::BufWriter::new) {
                if let Err(w) = serde_json::ser::to_writer_pretty(&mut wtr, &result.reads) {
                    debug!("{:?}", w);
                }
            }
            result
        };
        // Output gfa
        let filename = format!("{}/scaffolds.gfa", output_dir);
        let mut wtr = std::fs::File::create(&filename).map(std::io::BufWriter::new)?;
        writeln!(&mut wtr, "{}", result.gfa)?;
        result
            .assignments
            .into_iter()
            .filter_map(|(id, asn)| asn.map(|x| (id, x)))
            .collect()
    } else {
        use last_decompose::find_breakpoint::ReadClassify;
        encoded_reads
            .iter()
            .filter_map(|r| {
                initial_clusters
                    .iter()
                    .filter(|cl| cl.has(r.id()))
                    .nth(0)
                    .map(|cl| (r.id().to_string(), cl.id as u8))
            })
            .collect()
    };
    let mut decomposed: HashMap<u8, Vec<&fasta::Record>> = HashMap::new();
    let unassigned = results.values().copied().max().unwrap_or(0) + 1;
    for read in &reads {
        if let Some(cluster) = results.get(read.id()) {
            let cls = decomposed.entry(*cluster).or_insert(vec![]);
            cls.push(read);
        } else {
            let cls = decomposed.entry(unassigned).or_insert(vec![]);
            cls.push(read);
        }
    }
    let readlist = format!("{}/readlist.tsv", output_dir);
    let mut readlist = BufWriter::new(std::fs::File::create(readlist)?);
    let decomposed: HashMap<u8, Vec<_>> = decomposed
        .into_iter()
        //.filter(|(_, rs)| rs.len() > last_decompose::find_breakpoint::COVERAGE_THR)
        .collect();
    {
        for (&cluster_id, reads) in decomposed.iter() {
            let outpath = format!("{}/{}.fasta", output_dir, cluster_id);
            let wtr = match std::fs::File::create(&outpath) {
                Ok(res) => res,
                Err(why) => {
                    error!("{:?},{}", why, outpath);
                    continue;
                }
            };
            let mut wtr = fasta::Writer::new(wtr);
            for read in reads {
                let line = match read.desc() {
                    Some(d) => format!("{}\t{}\t{}", cluster_id, read.id(), d),
                    None => format!("{}\t{}\tNoDesc", cluster_id, read.id()),
                };
                writeln!(&mut readlist, "{}", line)?;
                wtr.write_record(read)?;
            }
        }
    }
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let dir = format!("{}/viewer", output_dir);
    let file = format!("{}/data.json", dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let res = dump_viewer(&results, &encoded_reads, &initial_clusters, &contigs)?;
    writeln!(&mut writer, "{}", res)?;
    Ok(())
}

fn create_viewer(matches: &clap::ArgMatches) -> std::io::Result<()> {
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
    let contigs = matches
        .value_of("contigs")
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
    let read_aln = matches
        .value_of("read_alignments")
        .map(|file| match last_tiling::parse_tab_file(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let contig_aln = matches
        .value_of("contig_alignments")
        .map(|file| match last_tiling::parse_tab_file(file) {
            Ok(res) => res,
            Err(why) => panic!("{}:{}", why, file),
        })
        .unwrap();
    let parse_line = |line: String| -> Option<(String, u8)> {
        let mut line = line.split("\t");
        let assign: u8 = line.next()?.parse().ok()?;
        let id = line.next()?.to_string();
        Some((id, assign))
    };
    let assignments: HashMap<String, u8> = matches
        .value_of("assignments")
        .map(|file| match std::fs::File::open(file) {
            Ok(res) => BufReader::new(res),
            Err(why) => panic!("{}:{}", why, file),
        })
        .map(|reader| {
            reader
                .lines()
                .filter_map(|e| e.ok())
                .filter_map(parse_line)
                .collect()
        })
        .unwrap();
    let thr: usize = matches
        .value_of("min_align_length")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let output_dir = matches.value_of("output_dir").unwrap();
    let dir = output_dir.to_string();
    if let Err(why) = std::fs::create_dir_all(&dir) {
        error!("Error Occured while outputing reads.");
        error!("{:?}", why);
        error!("This program did not work successfully.");
        error!("Shutting down...");
        std::process::exit(1);
    }
    use last_decompose::annotate_contigs_to_reference::annotate_aln_contigs_to_ref;
    let contigs_to_ref = annotate_aln_contigs_to_ref(&reference, &contig_aln, thr);
    let file = format!("{}/contig_alns.json", output_dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let contigs_to_ref = serde_json::ser::to_string(&contigs_to_ref)?;
    writeln!(&mut writer, "{}", contigs_to_ref)?;
    let contigs = last_tiling::Contigs::new(contigs);
    let reads = last_tiling::encoding(&reads, &contigs, &read_aln);
    use last_decompose::d3_data::convert_result_to_d3_data;
    let result_summary = convert_result_to_d3_data(&contigs, &reads, &assignments);
    let file = format!("{}/read_data.json", output_dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let result_summary = serde_json::ser::to_string(&result_summary)?;
    writeln!(&mut writer, "{}", result_summary)?;
    let file = format!("{}/linear.html", dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    writeln!(&mut writer, "{}", mito_assembler::template::TEMPLATE_LINEAR)?;
    Ok(())
}

fn main() -> std::io::Result<()> {
    let matches = App::new("MMMM")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Softwares to Decompose long reads.")
        .subcommand(subcommand_decompose())
        .subcommand(subcommand_create_viewer())
        .get_matches();
    match matches.subcommand() {
        ("decompose", Some(sub_m)) => decompose(sub_m),
        ("create_viewer", Some(sub_m)) => create_viewer(sub_m),
        _ => Ok(()),
    }
}
