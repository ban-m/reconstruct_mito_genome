extern crate bio;
extern crate bio_utils;
extern crate rusty_sandbox;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let sam_records = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    let fastq_records = rusty_sandbox::open_fastq_into_hashmap(&args[2])?;
    eprintln!("Opened dataset");
    let output = |target: String| -> std::io::Result<()> {
        let pileup = get_pileup_of(&sam_records[&target], &fastq_records, &target);
        eprintln!("Finish pileup. Pileup is {} height.", pileup.len());
        let mut file = loop {
            print!("Enter file to write:");
            std::io::stdout().flush()?;
            let file = get_input()?;
            match File::create(&Path::new(&file)).map(BufWriter::new) {
                Ok(res) => break res,
                Err(why) => eprintln!("Error:{:?}", why),
            };
        };
        let window = 200;
        let len = pileup[0].len();
        let mut start = 0;
        while start < len {
            let end = (start + window).min(len);
            for line in &pileup {
                writeln!(&mut file, "{}", String::from_utf8_lossy(&line[start..end]))?;
            }
            writeln!(&mut file)?;
            start += window;
        }
        // for line in pileup {
        //     writeln!(&mut file, "{}", String::from_utf8_lossy(&line))?;
        // }
        file.flush()
    };
    loop {
        print!("Enter read ID(or \"q\" to exit.):");
        std::io::stdout().flush()?;
        let target = get_input()?;
        if sam_records.contains_key(&target) {
            if let Err(why) = output(target) {
                eprintln!("{:?}", why);
                break;
            }
        } else if target == "q" {
            println!("Quitting...");
            break;
        } else if target == "LONGEST" {
            println!("Searching the longest read...");
            let target = fastq_records
                .iter()
                .filter(|(name, _)| sam_records.contains_key(name.as_str()))
                .max_by_key(|(_, seq)| seq.seq().len())
                .unwrap()
                .0
                .clone();
            println!("target:{}, {}", target, sam_records.contains_key(&target));
            if let Err(why) = output(target) {
                eprintln!("{:?}", why);
                break;
            }
        } else {
            println!("Invalid Read ID");
            println!("The IDs below are possible candidates");
            for key in sam_records.keys().take(10) {
                println!("{}", key);
            }
        }
    }
    Ok(())
}

fn get_input() -> std::io::Result<String> {
    let mut input = String::new();
    std::io::stdin().read_line(&mut input)?;
    Ok(input.trim().to_string())
}

use bio::io::fastq::Record;
use bio_utils::sam;
fn get_pileup_of(
    sams: &[sam::Sam],
    records: &HashMap<String, Record>,
    target: &str,
) -> Vec<Vec<u8>> {
    let coordinate: Vec<_> = records[target].seq().iter().copied().collect();
    let rlen = coordinate.len();
    let mut res:Vec<_> = sams.iter()
        .map(|sam| {
            let record = &records[sam.q_name()];
            let (start, end) = sam.mapped_region();
            let qlen = end - start;
            let pos = sam.pos();
            let mut rip = vec![b'_'; pos];
            rip.extend(&record.seq()[start..end]);
            if rlen > pos + qlen {
                rip.extend(vec![b'_'; rlen - pos - qlen])
            }
            rip
        })
        .chain(vec![coordinate])
        .collect();
    res.sort_by_key(|e| -(e.iter().take_while(|&&e| e == b'_').count() as isize));
    res
}
