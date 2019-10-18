extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate serde;
extern crate serde_json;
use env_logger::Env;
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
    // let fasta: Vec<_> = fasta
    //     .into_iter()
    //     .filter(|e| e.id() == "m54113_160913_184949/48432029/0_6399")
    //     .collect();
    use std::collections::HashMap;
    let readlen: HashMap<_, _> = fasta
        .iter()
        .map(|e| (e.id().to_string(), e.seq().len()))
        .collect();
    debug!("Read num\t{}", fasta.len());
    let repeats = last_tiling::repeat::open(&args[4])?;
    info!("Repeats:{:?}", repeats.len());
    let encoded_reads = last_tiling::encoding(&fasta, &contigs, &alignments);
    debug!("Encoded:\t{}", encoded_reads.len());
    for read in encoded_reads.iter() {
        let len = read.seq().iter().map(|e| e.len()).sum::<usize>();
        let raw_len = readlen[read.id()];
        if len != raw_len {
            eprintln!(
                "WARNING: READLEN DISCONCORDANCE {}\t{}\t{}",
                read.id(),
                len,
                raw_len
            );
        }
        // println!("{}", read);
    }
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
