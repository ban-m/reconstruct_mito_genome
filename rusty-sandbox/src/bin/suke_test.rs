extern crate bio;
extern crate rayon;
extern crate suke;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::time;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let fastq: Vec<_> = bio::io::fastq::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let k: usize = args[2].parse().unwrap();
    let n: usize = args[3].parse().unwrap();
    let reads: Vec<&[u8]> = fastq.iter().map(|e| e.seq()).collect();
    eprintln!("Building database...");
    let start = time::Instant::now();
    let multdb = suke::multihash::MultiHashDB::new(&reads, n, k);
    let singledb = suke::singlehash::SingleHashDB::new(&reads, n, k);
    let elapsed = time::Instant::now() - start;
    eprintln!("Database built!: {:?}", elapsed);
    let mut wtr = BufWriter::new(std::io::stdout());
    writeln!(&mut wtr, "ID1\tID2\tSingle\tMult")?;
    let start = time::Instant::now();
    for batch in fastq[..fastq.len()/2].chunks(10000) {
        let res:Vec<_> = batch
            .par_iter()
            .map(|x| {
                let maps:Vec<_> = multdb
                    .map(x.seq())
                    .into_iter()
                    .zip(singledb.map(x.seq()).into_iter())
                    .filter_map(|((y_idx, mscore), (y_idx2, sscore))| {
                        if mscore < 0.002 && sscore < 0.002 {
                            None
                        } else {
                            assert_eq!(y_idx, y_idx2);
                            Some((y_idx, mscore, sscore))
                        }
                    })
                    .collect::<Vec<_>>();
                (x, maps)
            })
            .collect();
        for (x, maps) in res {
            for (yidx, sim_single, sim_multi) in maps {
                writeln!(
                    &mut wtr,
                    "{}\t{}\t{}\t{}",
                    x.id(),
                    fastq[yidx].id(),
                    sim_single,
                    sim_multi
                )?;
            }
        }
    }
    eprintln!("Elapsed:{:?} ", time::Instant::now() - start);
    Ok(())
}
