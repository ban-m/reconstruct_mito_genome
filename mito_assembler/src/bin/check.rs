use std::collections::{HashMap, HashSet};
fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    // let alns = mito_assembler::last_alignment_train(&args[1], &args[2], 23).unwrap();
    // eprintln!("{}", alns.len());
    let chunked_reads: Vec<last_decompose::assemble::ChunkedRead> = std::fs::File::open(&args[1])
        .map(std::io::BufReader::new)
        .ok()
        .and_then(|e| serde_json::de::from_reader(e).ok())
        .unwrap();
    let labels: HashMap<_, HashSet<String>> =
        chunked_reads.iter().fold(HashMap::new(), |mut x, y| {
            if let Some(label) = y.label {
                x.entry(label).or_default().insert(y.id.clone());
            }
            x
        });
    let id2desc: HashMap<_, _> = chunked_reads
        .iter()
        .map(|read| (read.id.clone(), read.desc.clone()))
        .collect();
    let (k, thr) = (5, 15);
    let mut labels: Vec<_> = labels.into_iter().collect();
    labels.sort_by_key(|x| x.0);
    let result = last_decompose::assemble::assemble_reads(&chunked_reads, k, thr);
    for (id, asn) in result.0.iter() {
        match (id2desc[id].as_ref(), asn) {
            (Some(desc), Some(asn)) => println!("{}\t{}", asn, desc),
            (Some(desc), None) => println!("-\t{}", desc),
            (None, Some(asn)) => println!("{}\t-", asn),
            (None, None) => println!("-\t-"),
        }
    }
    let result: HashMap<_, HashSet<String>> =
        result.0.iter().fold(HashMap::new(), |mut x, (id, asn)| {
            if let Some(asn) = asn {
                x.entry(asn).or_default().insert(id.clone());
            }
            x
        });
    let mut result: Vec<_> = result.into_iter().collect();
    result.sort_by_key(|x| x.0);
    for (init_cl, members) in labels.iter() {
        eprintln!("{}\t{}", init_cl, members.len());
    }
}
