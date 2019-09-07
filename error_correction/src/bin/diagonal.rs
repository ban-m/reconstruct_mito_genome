extern crate bio;
extern crate edlib_sys;
extern crate libc;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rayon;
use rayon::prelude::*;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let result: Vec<_> = {
        debug!("start");
        let before: HashMap<_, _> = bio::io::fastq::Reader::from_file(&args[1])?
            .records()
            .filter_map(|e| e.ok())
            .filter_map(|rec| {
                let id: usize = rec.id().parse().ok()?;
                let seq = rec.seq().to_vec();
                Some((id, seq))
            })
            .collect();
        let corrected: Vec<(usize, _)> = bio::io::fasta::Reader::from_file(&args[2])?
            .records()
            .filter_map(|e| e.ok())
            .filter_map(|rec| {
                let id: usize = rec.id().split('_').nth(0)?.parse().ok()?;
                let seq = rec.seq().to_vec();
                Some((id, seq))
            })
            .collect();
        debug!("open records.");
        corrected
            .into_par_iter()
            .filter_map(|(id, query)| {
                before
                    .get(&id)
                    .map(|target| (id, edlib_sys::align(&query, &target)))
            })
            .collect()
    };
    debug!("Output");
    use std::io::{BufWriter, Write};
    let out = std::io::stdout();
    let mut out = BufWriter::new(out.lock());
    for (id, cigar) in result {
        out.write_all(id.to_string().as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(cigar.as_bytes())?;
        out.write_all(b"\n")?;
    }
    Ok(())
}
