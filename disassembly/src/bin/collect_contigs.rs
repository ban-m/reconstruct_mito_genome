extern crate bio_utils;
use bio_utils::fasta;
use std::fs;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::path::PathBuf;
    let contigs: Vec<(usize, PathBuf)> = fs::read_dir(&args[1])?
        .filter_map(|e| e.ok())
        .filter_map(|entry| {
            let index = entry.path();
            let index = index.file_name()?.to_str()?.to_string();
            let mut contig = std::path::PathBuf::new();
            contig.push(&args[1]);
            contig.push(&index);
            contig.push("scaffolds.fasta");
            let index: usize = index.parse().ok()?;
            Some((index, contig))
        })
        .collect();
    let contigs: Vec<fasta::Record> = contigs
        .into_iter()
        .filter_map(|(index, path)| {
            let contig: Vec<_> = fasta::Reader::from_file(&path)
                .ok()?
                .records()
                .filter_map(|e| e.ok())
                .enumerate()
                .map(|(num, record)| {
                    let id = format!("{}_{}", index, num);
                    fasta::Record::with_data(&id, &None, record.seq())
                })
                .collect();
            let mut wtr_path = PathBuf::new();
            wtr_path.push(&args[2]);
            wtr_path.push(&format!("{}.contigs.fasta", index));
            let wtr = fs::File::create(wtr_path).ok()?;
            let mut wtr = fasta::Writer::new(wtr);
            for r in contig.iter() {
                wtr.write_record(r).ok()?
            }
            Some(contig)
        })
        .flat_map(|c| c)
        .collect();
    let wtr = format!("{}/multipartite.fasta", &args[2]);
    let wtr = fs::File::create(wtr)?;
    let mut wtr = fasta::Writer::new(wtr);
    for c in contigs {
        wtr.write_record(&c)?;
    }
    Ok(())
}
