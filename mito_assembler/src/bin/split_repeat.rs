extern crate bio_utils;
extern crate last_tiling;
use bio_utils::fasta;
const THR: usize = 1_000;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut contigs: Vec<_> = fasta::parse_into_vec(&args[1])?;
    let mut alns: Vec<_> = last_tiling::parse_tab_file(&args[2])?
        .into_iter()
        .filter(|e| e.seq1_matchlen() > THR && e.seq2_matchlen() > THR)
        .collect();
    let mut result = vec![];
    while !alns.is_empty() {
        let count = count_repetitiveness(&contigs, &alns);
        let (max_contig, max_position) = find_maximum(&count);
        let (mut start, mut end) = (max_position, max_position);
        let c = &count[max_contig];
        while 0 < c[start] && 0 < start {
            start -= 1;
        }
        while 0 < c[end] && end < c.len() {
            end += 1;
        }
        let contig = contigs
            .iter_mut()
            .filter(|e| e.id() == max_contig)
            .nth(0)
            .unwrap();
        result.push(fasta::Record::with_data(
            &format!("rep{}", result.len()),
            &None,
            &contig.seq()[start..end],
        ));
        alns = alns
            .into_iter()
            .filter(|aln| {
                (aln.seq1_name() == max_contig
                    && start < aln.seq1_start_from_forward()
                    && aln.seq1_end_from_forward() < end)
                    || (aln.seq2_name() == max_contig
                        && start < aln.seq2_start_from_forward()
                        && aln.seq2_end_from_forward() < end)
            })
            .collect::<Vec<_>>();
        let seq: Vec<_> = contig.seq()[..start]
            .iter()
            .chain(contig.seq()[end..].iter())
            .copied()
            .collect();
        *contig = fasta::Record::with_data(contig.id(), contig.desc(), &seq);
    }
    let stdout = std::io::stdout();
    use std::io::BufWriter;
    let mut wtr = fasta::Writer::new(BufWriter::new(stdout.lock()));
    for record in contigs
        .into_iter()
        .enumerate()
        .map(|(idx, r)| fasta::Record::with_data(&format!("{}", idx), &None, r.seq()))
        .chain(result.into_iter())
    {
        wtr.write_record(&record)?;
    }
    Ok(())
}

fn find_maximum<'a>(count: &'a HashMap<String, Vec<usize>>) -> (&'a str, usize) {
    let (x, y) = count
        .iter()
        .filter_map(|(key, val)| {
            let max = val.iter().enumerate().max_by_key(|(_, x)| *x)?;
            Some((key, max.0))
        })
        .max_by_key(|(_, val)| *val)
        .unwrap();
    (x, y)
}

fn count_repetitiveness(
    contigs: &Vec<fasta::Record>,
    alns: &[last_tiling::LastTAB],
) -> HashMap<String, Vec<usize>> {
    let mut count: HashMap<_, _> = contigs
        .iter()
        .map(|e| (e.id().to_string(), vec![0; e.seq().len()]))
        .collect();
    for aln in alns.iter() {
        let c1 = count.get_mut(aln.seq1_name()).unwrap();
        for i in aln.seq1_start_from_forward()..aln.seq1_end_from_forward() {
            c1[i] += 1;
        }
        let c2 = count.get_mut(aln.seq2_name()).unwrap();
        for i in aln.seq2_start_from_forward()..aln.seq2_end_from_forward() {
            c2[i] += 1;
        }
    }
    count
}
