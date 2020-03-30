extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate mito_assembler;
extern crate serde_json;
use clap::{App, Arg};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
fn main() -> std::io::Result<()> {
    let matches = App::new("M4_CreateViewer")
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
                .help("Assignmnet of each read")
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
    let output_dir = matches.value_of("output_dir").unwrap();
    let dir = format!("{}", output_dir);
    if let Err(why) = std::fs::create_dir_all(&dir) {
        error!("Error Occured while outputing reads.");
        error!("{:?}", why);
        error!("This program did not work successfully.");
        error!("Shutting down...");
        std::process::exit(1);
    }
    use last_decompose::annotate_contigs_to_reference::annotate_aln_contigs_to_ref;
    let contigs_to_ref = annotate_aln_contigs_to_ref(&reference, &contig_aln);
    let file = format!("{}/contig_alns.json", output_dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let contigs_to_ref = serde_json::ser::to_string(&contigs_to_ref)?;
    writeln!(&mut writer, "{}", contigs_to_ref)?;

    let contigs = last_tiling::Contigs::new(contigs.clone());
    let reads = last_tiling::encoding(&reads, &contigs, &read_aln);
    use last_decompose::d3_data::convert_result_to_d3_data;
    let result_summary = convert_result_to_d3_data(&contigs, &reads, &assignments);
    let file = format!("{}/read_data.json", output_dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    let result_summary = serde_json::ser::to_string(&result_summary)?;
    writeln!(&mut writer, "{}", result_summary)?;
    Ok(())
}
