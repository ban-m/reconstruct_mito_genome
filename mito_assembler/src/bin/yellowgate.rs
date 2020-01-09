extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate mito_assembler;
use mito_assembler::dump_viewer;
use clap::{App, Arg};
use last_decompose::find_breakpoint::{ReadClassify};
use last_decompose::ERead;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let matches = App::new("YellowGate")
        .version("0.1")
        .author("Bansho MASUTANI")
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
        .get_matches();
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
    let config = last_decompose::error_profile::summarize_tab(&alignments, &reads, &reference);
    debug!("Profiled Error Rates:{}", config);
    let contigs = last_tiling::contig::Contigs::new(reference);
    let repeats = last_tiling::into_repeats(&self_aln, &contigs);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let encoded_reads: Vec<_> = encoded_reads
        .into_iter()
        .map(ERead::new_no_gapfill)
        .collect();
    let critical_regions = last_decompose::critical_regions(&encoded_reads, &contigs, &repeats);
    for c in &critical_regions {
        let counts = encoded_reads
            .iter()
            .filter(|read| c.along_with(read))
            .count();
        debug!("{:?} having {} reads.", c, counts);
    }
    use std::collections::HashMap;
    let cr = &critical_regions;
    let mock_ans: HashMap<String, u8> = HashMap::new();
    let results: HashMap<String, u8> =
        last_decompose::decompose(encoded_reads, &cr, &contigs, &repeats, config, &mock_ans)
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
    for (&cluster_id, reads) in decomposed.iter() {
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
    let res = dump_viewer(&results, &encoded_reads, &critical_regions, &contigs)?;
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
    writeln!(&mut writer, "{}", res)?;
    Ok(())
}
