extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use env_logger::Env;
use std::io::Write;
fn main() -> std::io::Result<()> {
    // env_logger::from_env(Env::default().default_filter_or("debug")).init();
    env_logger::from_env(Env::default().default_filter_or("warn")).init();
    let args: Vec<_> = std::env::args().collect();
    info!("Start");
    let alignments = last_tiling::parse_tab_file(&args[1])?;
    debug!("Alignments:{}", alignments.len());
    let contigs = last_tiling::contig::Contigs::from_file(&args[2])?;
    debug!("Contig files:\n{}", contigs);
    let fasta = bio_utils::fasta::parse_into_vec(&args[3])?;
    debug!("Read num\t{}", fasta.len());
    let encoded_reads = last_tiling::encoding(&fasta, &contigs, &alignments);
    debug!("Encoded:\t{}", encoded_reads.len());
    // let out = std::io::stdout();
    // let mut out = BufWriter::new(out.lock());
    // for read in &encoded_reads {
    //     writeln!(&mut out, "{}", read)?;
    // }
    eprintln!("Output dump");
    let mut wtr = std::fs::File::create("./data/contigs.json")?;
    wtr.write_all(serde_json::ser::to_string_pretty(&contigs)?.as_bytes())?;
    let mut wtr = std::fs::File::create("./data/reads.json")?;
    wtr.write_all(serde_json::ser::to_string(&encoded_reads)?.as_bytes())?;
    // for read in encoded_reads {
    //     println!("{:?}", read.id);
    //     for unit in &read.seq {
    //         match unit {
    //             last_tiling::unit::ChunkedUnit::En(encode) => {
    //                 let start = encode.unit as usize * last_tiling::UNIT_SIZE;
    //                 let end = (encode.unit + 1) as usize * last_tiling::UNIT_SIZE;
    //                 assert!(start < end,"{},{},{:?}",start,end,read);
    //                 let refr = if encode.is_forward() {
    //                     println!("[F]");
    //                     &contigs.get_by_id(encode.contig).unwrap()[start..end]
    //                 } else {
    //                     println!("[R]");
    //                     let len = contigs.get_by_id(encode.contig).unwrap().len();
    //                     let (start,end) = (len - end, len -start);
    //                     assert!(start <end,"{},{},{}",start,end,len);
    //                     &contigs.get_by_id_revcmp(encode.contig).unwrap()[start..end]
    //                 };
    //                 encode.view(refr);
    //             }
    //             last_tiling::unit::ChunkedUnit::Gap(gap) => {
    //                 println!("Gap:{:?}", gap);
    //             }
    //         }
    //     }
    // }
    // let contig = 0;
    // let unit_n = 20;
    // for read in encoded_reads {
    //     for unit in &read.seq {
    //         match unit {
    //             last_tiling::unit::ChunkedUnit::En(encode) => {
    //                 if encode.contig == contig && encode.unit == unit_n {
    //                     let start = encode.unit as usize * last_tiling::UNIT_SIZE;
    //                     let end = (encode.unit + 1) as usize * last_tiling::UNIT_SIZE;
    //                     let refr = if encode.is_forward() {
    //                         println!("[F]");
    //                         &contigs.get_by_id(encode.contig).unwrap()[start..end]
    //                     } else {
    //                         println!("[R]");
    //                         let len = contigs.get_by_id(encode.contig).unwrap().len();
    //                         let (start, end) = (len - end, len - start);
    //                         assert!(start < end, "{},{},{}", start, end, len);
    //                         &contigs.get_by_id_revcmp(encode.contig).unwrap()[start..end]
    //                     };
    //                     encode.view(refr);
    //                 }
    //             }
    //             _ => {}
    //         }
    //     }
    // }
    Ok(())
}
