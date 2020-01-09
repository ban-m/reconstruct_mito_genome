extern crate bio_utils;
extern crate last_decompose;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let reads = match bio_utils::fasta::parse_into_vec(&args[1]) {
        Ok(res) => res,
        Err(why) => panic!("{}:{}", why, &args[1]),
    };
    let reference = match bio_utils::fasta::parse_into_vec(&args[2]) {
        Ok(res) => res,
        Err(why) => panic!("{}:{}", why, &args[2]),
    };
    let alignments = match last_tiling::parse_tab_file(&args[3]) {
        Ok(res) => res,
        Err(why) => panic!("{}:{}", why, &args[3]),
    };
    let self_aln = match last_tiling::parse_tab_file(&args[4]) {
        Ok(res) => res,
        Err(why) => panic!("{}:{}", why, &args[4]),
    };
    let output_dir = &args[5];
    let answer: HashMap<String, u8> = reads
        .iter()
        .map(|e| {
            let id = e.id().to_string();
            let desc = match e.desc() {
                Some(res) => res,
                None => return (id, 0),
            };
            if desc.contains("master_circle") {
                (id, 0)
            } else {
                (id, 1)
            }
        })
        .collect();
    let config = last_decompose::error_profile::summarize_tab(&alignments, &reads, &reference);
    debug!("Profiled Error Rates:{}", config);
    let contigs = last_tiling::contig::Contigs::new(reference);
    let repeats = last_tiling::into_repeats(&self_aln, &contigs);
    let encoded_reads = last_tiling::encoding(&reads, &contigs, &alignments);
    use last_decompose::ERead;
    let encoded_reads: Vec<_> = encoded_reads
        .into_iter()
        .map(ERead::new_no_gapfill)
        .collect();
    let critical_regions = last_decompose::critical_regions(&encoded_reads, &contigs, &repeats);
    for c in &critical_regions {
        use last_decompose::find_breakpoint::ReadClassify;
        let counts = encoded_reads
            .iter()
            .filter(|read| c.along_with(read))
            .count();
        debug!("{:?} having {} reads.", c, counts);
    }
    use std::collections::HashMap;
    let cr = &critical_regions;
    let results: HashMap<String, u8> =
        last_decompose::decompose(encoded_reads, cr, &contigs, &repeats, config, &answer)
            .into_iter()
            .collect();
    use bio_utils::fasta;
    let mut decomposed: HashMap<u8, Vec<fasta::Record>> = HashMap::new();
    for read in reads {
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
    for (cluster_id, reads) in decomposed {
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
            wtr.write_record(&read)?;
        }
    }
    Ok(())
}
