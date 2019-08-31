extern crate bio;
extern crate libc;
extern crate rand;
extern crate rayon;
extern crate select_mitochondrial_genome;
use bio::io::fastq::Record;
use rand::{seq::SliceRandom, thread_rng};
use rayon::prelude::*;
use select_mitochondrial_genome::ReadSummary;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
fn open_records(file: &str) -> std::io::Result<Vec<Record>> {
    let mut records: Vec<_> = bio::io::fastq::Reader::from_file(&Path::new(&file))?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let mut rng = thread_rng();
    records.shuffle(&mut rng);
    Ok(records)
}

fn create_fifo(name: String) -> String {
    let cfile = std::ffi::CString::new(name.clone()).unwrap();
    unsafe {
        libc::mkfifo(cfile.as_ptr(), 0o644);
    }
    name
}

fn write_retry<W: std::io::Write>(record: &Record, wtr: &mut BufWriter<W>) {
    let mut output = vec![b'>'];
    output.extend(record.id().as_bytes());
    output.push(b'\n');
    output.extend(record.seq());
    output.push(b'\n');
    let mut is_cond = 0;
    while let Err(_) = wtr.write_all(&output) {
        is_cond += 1;
    }
    if is_cond != 0 {
        eprintln!("Write {} sucess in {} trial", record.id(), is_cond);
    }
}

// return paf.
fn all_vs_all_alignments(records: &[Record], prefix: &str) -> std::io::Result<String> {
    let fifo1 = create_fifo(format!("{}-1", prefix));
    let fifo2 = create_fifo(format!("{}-2", prefix));
    let child = std::process::Command::new("minimap2")
        .args(&["-x", "ava-pb"])
        .args(&["-t", "1"])
        .arg(&fifo1)
        .arg(&fifo2)
        .stdout(std::process::Stdio::piped())
        .spawn()?;
    let mut rm = std::process::Command::new("rm");
    rm.arg(&fifo1).arg(&fifo2);
    let open = |name| -> std::io::Result<_> {
        Ok(BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .read(false)
                .create(false)
                .open(name)?,
        ))
    };
    if let Ok(mut input1) = open(fifo1) {
        for record in records {
            write_retry(record, &mut input1);
        }
    }
    if let Ok(mut input2) = open(fifo2) {
        for record in records {
            write_retry(record, &mut input2);
        }
    }
    let out = child.wait_with_output().expect("failed to get out");
    rm.output().expect("");
    Ok(String::from_utf8_lossy(&out.stdout).to_string())
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let records = open_records(&args[1])?;
    let coverage: usize = args[2].parse().unwrap();
    let output_prefix = &args[3];
    let chunksize = records.len() / coverage;
    let chunks: Vec<_> = (0..coverage)
        .map(|i| {
            let start = chunksize * i;
            let end = if i == coverage - 1 {
                records.len()
            } else {
                chunksize * (i + 1)
            };
            (start, end)
        })
        .collect();
    let (summaries, pafs): (Vec<(String, ReadSummary)>, Vec<String>) = chunks
        .into_par_iter()
        .enumerate()
        .filter_map(|(idx,(s, e))| {
            let chunk = &records[s..e];
            let prefix = format!("{}.{}.temp",output_prefix,idx);
            let paf = all_vs_all_alignments(chunk, &prefix).ok()?;
            let summaries: Vec<_> = select_mitochondrial_genome::string_to_intervals(&paf)
                .into_iter()
                .map(|(id, interval)| (id, ReadSummary::from_interval(interval)))
                .collect();
            Some((summaries, paf))
        })
        .fold(
            || (vec![], vec![]),
            |(mut sums, mut pafs), (sum, paf)| {
                sums.extend(sum);
                pafs.push(paf);
                (sums, pafs)
            },
        )
        .reduce(
            || (vec![], vec![]),
            |(mut sums, mut pafs), (x, y)| {
                sums.extend(x);
                pafs.extend(y);
                (sums, pafs)
            },
        );
    let summary_wtr = format!("{}.tsv", output_prefix);
    let paf_wtr = format!("{}.paf", output_prefix);
    let mut summary_wtr = BufWriter::new(File::create(&Path::new(&summary_wtr))?);
    let mut paf_wtr = BufWriter::new(File::create(&Path::new(&paf_wtr))?);
    for (id, s) in summaries {
        writeln!(
            &mut summary_wtr,
            "{}\t{}\t{}\t{}",
            id, s.mean, s.sd, s.length
        )?;
    }
    for paf in pafs {
        write!(&mut paf_wtr, "{}", paf)?;
    }
    Ok(())
}
