extern crate bio_utils;
extern crate env_logger;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate serde;
extern crate serde_json;
use last_decompose::ERead;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let alignments = last_tiling::parse_tab_file(&args[2])?;
    let contigs = last_tiling::Contigs::from_file(&args[3])?;
    let reads = last_tiling::encoding(&reads, &contigs, &alignments);
    let reads: Vec<_> = reads.into_iter().map(ERead::new_no_gapfill).collect();
    let res = last_decompose::find_breakpoint::confluent_position(&alignments, &contigs, 150);
    for e in res {
        debug!("{:?}", e);
    }
    let res = last_decompose::find_breakpoint::contigpair_position(&reads, &contigs);
    for e in res {
        debug!("{:?}", e);
    }
    return Ok(());
    // let contigs = bio_utils::fasta::parse_into_vec(&args[3])?;
    // let repeats = last_tiling::repeat::open(&args[4])?;
    // debug!("Start");
    // let contigs = last_tiling::Contigs::new(contigs);
    // let reads: Vec<_> = last_tiling::encoding(&reads, &contigs, &alignments)
    //     .into_iter()
    //     //.inspect(|e| debug!("{}", e))
    //     .map(ERead::new)
    //     //.inspect(|e| debug!("{}", e))
    //     .collect();
    // debug!("Start");
    // let result = last_decompose::critical_regions(&reads, &contigs, &repeats);
    // use last_decompose::find_breakpoint::ReadClassify;
    // for cr in &result {
    //     let num = reads.iter().filter(|r| cr.along_with(r)).count();
    //     debug!("{:?} has {} reads.", cr, num);
    // }
    // let wtr = std::io::stdout();
    // let mut wtr = std::io::BufWriter::new(wtr.lock());
    // serde_json::ser::to_writer_pretty(&mut wtr, &result).unwrap();
    // debug!("Enumerating all co-occurence of Critical Regions");
    // use last_decompose::CriticalRegion;
    // use std::collections::HashMap;
    // let result: HashMap<Vec<&CriticalRegion>, u32> = reads
    //     .iter()
    //     .filter_map(|read| {
    //         let mut along_crs: Vec<&CriticalRegion> =
    //             result.iter().filter(|cr| cr.along_with(read)).collect();
    //         along_crs.sort();
    //         if along_crs.len() > 1 {
    //             Some(along_crs)
    //         } else {
    //             None
    //         }
    //     })
    //     .fold(HashMap::new(), |mut current, x| {
    //         current.entry(x).and_modify(|e| *e += 1).or_insert(1);
    //         current
    //     });
    // info!("Dump Critical regions which co-occured in reads with counts.");
    // for (crs, count) in result.into_iter().filter(|&(_, count)| count > 5) {
    //     info!("The CR below appeared in {} times in the dataset.", count);
    //     for cr in crs {
    //         eprintln!("{}  ->", cr);
    //     }
    //     info!("");
    // }
    // Ok(())
}
