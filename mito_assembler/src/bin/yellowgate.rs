extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
use clap::{App, Arg};
fn main() -> std::io::Result<()> {
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
                .long("contig")
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
    let contigs = last_tiling::contig::Contigs::new(reference);
    let repeats = last_tiling::into_repeats(&self_aln, &contigs);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    use last_decompose::ERead;
    let encoded_reads: Vec<_> = encoded_reads.into_iter().map(ERead::new).collect();
    let critical_regions = last_decompose::critical_regions(&encoded_reads, &contigs, &repeats);
    use std::collections::HashMap;
    let results: HashMap<String, u8> =
        last_decompose::decompose(encoded_reads, &critical_regions, &contigs, &repeats, config)
            .into_iter()
            .collect();
    use bio_utils::fasta;
    let mut decomposed: HashMap<u8, Vec<fasta::Record>> = HashMap::new();
    for read in reads {
        if let Some(cluster) = results.get(read.id()) {
            let cls = decomposed.entry(*cluster).or_insert(vec![]);
            cls.push(read);
        } else {
            eprint!("Read:{} is not in the results. ", read.id());
            eprintln!("Maybe it is a junk read. Go to the next read.");
        }
    }
    if let Err(why) = std::fs::create_dir_all(output_dir) {
        eprintln!("Error Occured while outputing reads.");
        eprintln!("{:?}", why);
        eprintln!("This program did not work successfully.");
        eprintln!("Shutting down...");
        std::process::exit(1);
    }
    for (cluster_id, reads) in decomposed {
        let outpath = format!("{}/{}.fasta", output_dir, cluster_id);
        let wtr = match std::fs::File::create(&outpath) {
            Ok(res) => res,
            Err(why) => {
                eprintln!("Error Occured while creating a file:{}", outpath);
                eprintln!("{:?}", why);
                eprintln!("This program did not work successfully.");
                eprintln!("Shutting down...");
                std::process::exit(1);
            }
        };
        let mut wtr = fasta::Writer::new(wtr);
        for read in reads {
            wtr.write_record(&read).unwrap()
        }
    }
    Ok(())
}
