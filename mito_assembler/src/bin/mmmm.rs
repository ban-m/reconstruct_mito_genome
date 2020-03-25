extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
extern crate serde_json;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate mito_assembler;
use clap::{App, Arg};
use last_decompose::ERead;
use mito_assembler::dump_viewer;
fn main() -> std::io::Result<()> {
    let matches = App::new("MMMM")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Decomposing long reads.")
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
        .get_matches();
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("MMMM started. Debug mode.");
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
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    debug!("Profiled Error Rates:{}", config);
    let contigs = last_tiling::contig::Contigs::new(reference);
    let repeats = last_tiling::into_repeats(&self_aln, &contigs);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let start_stop = last_tiling::get_start_stop(&encoded_reads, &contigs);
    let encoded_reads: Vec<_> = encoded_reads
        .into_iter()
        .map(ERead::new_no_gapfill)
        .collect();
    let initial_clusters =
        last_decompose::initial_clusters(&encoded_reads, &contigs, &repeats, &alignments);
    debug!("Initial clusters constructed");
    for c in &initial_clusters {
        let counts = c.ids().len();
        debug!("{:?} having {} reads.", c, counts);
    }
    use std::collections::HashMap;
    let mock_ans: HashMap<String, u8> = HashMap::new();
    let cl = cluster_num;
    let results: HashMap<String, u8> = last_decompose::decompose(
        encoded_reads,
        &initial_clusters,
        &contigs,
        &repeats,
        &config,
        &mock_ans,
        cl,
    )
    .into_iter()
    .collect();
    use bio_utils::fasta;
    let mut decomposed: HashMap<u8, Vec<&fasta::Record>> = HashMap::new();
    for read in &reads {
        if let Some(cluster) = results.get(read.id()) {
            let cls = decomposed.entry(*cluster).or_insert(vec![]);
            cls.push(read);
        } else {
            warn!("Read:{} is not in the results. ", read.id());
            warn!("Maybe it is a junk read. Go to the next read.");
        }
    }
    if let Err(why) = std::fs::create_dir_all(output_dir) {
        error!("Error Occured while outputing reads.");
        error!("{:?}", why);
        error!("This program did not work successfully.");
        error!("Shutting down...");
        std::process::exit(1);
    }
    use std::io::{BufWriter, Write};
    let readlist = format!("{}/readlist.tsv", output_dir);
    let mut readlist = BufWriter::new(std::fs::File::create(readlist)?);
    for (&cluster_id, reads) in decomposed
        .iter()
        .filter(|&(_, rs)| rs.len() > last_decompose::find_breakpoint::COVERAGE_THR)
    {
        let outpath = format!("{}/{}.fasta", output_dir, cluster_id);
        let wtr = match std::fs::File::create(&outpath) {
            Ok(res) => res,
            Err(why) => {
                error!("Error Occured while creating a file:{}", outpath);
                error!("{:?}", why);
                error!("This program did not work successfully.");
                error!("Shutting down...");
                std::process::exit(1);
            }
        };
        let mut wtr = fasta::Writer::new(wtr);
        for read in reads {
            writeln!(&mut readlist, "{}\t{}", cluster_id, read.id())?;
            wtr.write_record(read)?;
        }
    }
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let file = format!("{}/start_stop.tsv", output_dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    for (contig, stst) in start_stop {
        for (pos, count) in stst {
            writeln!(&mut writer, "{}\t{}\t{}", contig, pos, count)?;
        }
    }
    let dir = format!("{}/viwer", output_dir);
    if let Err(why) = std::fs::create_dir_all(&dir) {
        error!("Error Occured while outputing reads.");
        error!("{:?}", why);
        error!("This program did not work successfully.");
        error!("Shutting down...");
        std::process::exit(1);
    }
    let file = format!("{}/data.json", dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let res = dump_viewer(&results, &encoded_reads, &initial_clusters, &contigs)?;
    writeln!(&mut writer, "{}", res)?;
    let file = format!("{}/repeats.json", dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let repeats = serde_json::ser::to_string(&repeats).unwrap();
    writeln!(&mut writer, "{}", repeats)?;
    Ok(())
}
