use bio_utils::fasta;
use bio_utils::maf;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let contigs: Vec<_> = fasta::parse_into_vec(&args[1])?;
    let mut coverages: HashMap<_, Vec<u32>> = contigs
        .iter()
        .map(|e| (e.id().to_string(), vec![0; e.seq().len()]))
        .collect();
    let alignments: Vec<_> = maf::Reader::from_file(&args[2])?
        .records()
        .filter_map(|e| e.ok())
        .collect();
    for aln in alignments.into_iter() {
        let start = aln.sequence()[0].start() as usize;
        let coverage = coverages.get_mut(aln.sequence()[0].name()).unwrap();
        let mut pos = start;
        for (&rfr, &qry) in aln.sequence()[0]
            .text()
            .iter()
            .zip(aln.sequence()[1].text().iter())
        {
            pos += (rfr != b'-') as usize;
            coverage[pos - 1] += (qry != b'-') as u32;
        }
    }
    let out = std::io::stdout();
    let mut wtr = fasta::Writer::new(std::io::BufWriter::new(out.lock()));
    for record in contigs {
        let coverage = coverages.get(record.id()).unwrap();
        let len = coverage.len() as u32;
        let mean = coverage.iter().sum::<u32>() / len;
        let sqave = coverage.iter().map(|&c| (c - mean).pow(2u32)).sum::<u32>() / len;
        let sd = ((sqave - mean * mean) as f64).sqrt().floor() as u32;
        let thr = (mean.max(3 * sd) - 3 * sd).max(5);
        eprintln!("MEAN\tSD\tTHR:{}\t{}\t{}", mean, sd, thr);
        let (start, stop) = get_range(coverage, thr);
        let seq: Vec<_> = record.seq()[start..stop].to_vec();
        let desc = record.desc().cloned();
        let record = fasta::Record::with_data(record.id(), &desc, &seq);
        wtr.write_record(&record)?;
    }
    Ok(())
}

fn get_range(coverage: &[u32], thr: u32) -> (usize, usize) {
    // Start position.
    let mut cum_score = vec![];
    let len = coverage.len() as f64;
    let mut score = 0.;
    for &c in coverage.iter() {
        score += if c < thr { len.recip() } else { -len.recip() };
        cum_score.push(score);
    }
    let (argmax, max) = cum_score
        .iter()
        .enumerate()
        .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
        .unwrap();
    eprintln!("START:{}\t{}", max, argmax);
    let start = if *max < 0. { 0 } else { argmax };
    // End position.
    score = 0.;
    for (idx, &c) in coverage.iter().enumerate().rev() {
        score += if c < thr { len.recip() } else { -len.recip() };
        cum_score[idx] = score;
    }
    let (argmax, max) = cum_score
        .iter()
        .enumerate()
        .skip(start)
        .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
        .unwrap();
    eprintln!("END:{}\t{}", max, argmax);
    let end = if *max < 0. { coverage.len() } else { argmax };
    (start, end)
}
