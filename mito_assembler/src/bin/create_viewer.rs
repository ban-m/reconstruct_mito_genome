extern crate bio_utils;
extern crate clap;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate mito_assembler;
extern crate serde_json;
use last_decompose::ERead;
use mito_assembler::dump_viewer;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let reference = bio_utils::fasta::parse_into_vec(&args[2])?;
    let alignments = last_tiling::parse_tab_file(&args[3])?;
    let self_aln = last_tiling::parse_tab_file(&args[4])?;
    let output_dir = &args[5];
    use std::io::{BufRead, BufReader, BufWriter, Write};
    let results: HashMap<String, u8> = BufReader::new(std::fs::File::open(&args[6])?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let e: Vec<_> = e.split('\t').collect::<Vec<_>>();
            let cluster: u8 = e[0].parse().ok()?;
            let id = e[1].to_string();
            Some((id, cluster))
        })
        .collect();
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
    let file = format!("{}/repeats.json", dir);
    let mut writer = BufWriter::new(std::fs::File::create(&file)?);
    serde_json::ser::to_writer_pretty(&mut writer, &repeats).unwrap();
    Ok(())
}
