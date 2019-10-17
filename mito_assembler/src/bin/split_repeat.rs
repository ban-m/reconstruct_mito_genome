extern crate bio_utils;
extern crate last_tiling;
#[macro_use]
extern crate log;
extern crate env_logger;
use bio_utils::fasta;
const THR: usize = 1_000;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    use env_logger::Env;
    env_logger::from_env(Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let contigs: Vec<_> = fasta::parse_into_vec(&args[1])?;
    let mut alns = {
        let alns: Vec<_> = last_tiling::parse_tab_file(&args[2])?
            .into_iter()
            .filter(|e| e.seq1_matchlen() > THR && e.seq2_matchlen() > THR)
            .filter(|e| {
                e.seq1_len() != e.seq1_matchlen()
                    && e.seq2_len() != e.seq2_matchlen()
                    && e.seq1_name() != e.seq2_name()
            })
            .collect();
        let mut remove_dedup: Vec<last_tiling::LastTAB> = vec![];
        for aln in alns {
            if remove_dedup.iter().all(|a| &aln != a) {
                remove_dedup.push(aln);
            }
        }
        let mut remove_contained = vec![];
        debug!("{}", remove_dedup.len());
        for aln in &remove_dedup {
            if remove_dedup.iter().all(|a| {
                if aln == a {
                    true
                } else {
                    let (s1_name, s2_name) = (aln.seq1_name(), aln.seq2_name());
                    let (c1_name, c2_name) = (a.seq1_name(), a.seq2_name());
                    let (s1_s, s1_e) = (aln.seq1_start_from_forward(), aln.seq1_end_from_forward());
                    let (s2_s, s2_e) = (aln.seq2_start_from_forward(), aln.seq2_end_from_forward());
                    let (c1_s, c1_e) = (a.seq1_start_from_forward(), a.seq1_end_from_forward());
                    let (c2_s, c2_e) = (a.seq2_start_from_forward(), a.seq2_end_from_forward());
                    // 1 contains 2.
                    let is_cont = |name1, name2, start1, start2, end1, end2| -> bool {
                        name1 == name2 && start1 <= start2 && end2 <= end1
                    };
                    !is_cont(c1_name, s1_name, c1_s, s1_s, c1_e, s1_e)
                        && !is_cont(c1_name, s2_name, c1_s, s2_s, c1_e, s2_e)
                        && !is_cont(c2_name, s1_name, c2_s, s1_s, c2_e, s1_e)
                        && !is_cont(c2_name, s2_name, c2_s, s2_s, c2_e, s2_e)
                }
            }) {
                remove_contained.push(aln.clone());
            }
        }
        remove_contained
    };
    let mut result = vec![];
    for c in &contigs {
        debug!("{}:{}", c.id(), c.seq().len());
    }
    // If mask['contig'][i] = false, it would be trimmed as repetitive region.
    let mask: HashMap<_, _> = count_repetitiveness(&contigs, &alns)
        .into_iter()
        .map(|(key, val)| (key, val.into_iter().map(|e| e == 0).collect::<Vec<_>>()))
        .collect();
    while !alns.is_empty() {
        debug!("remaining repeats:{}", alns.len());
        for aln in &alns {
            debug!(
                "{}:{}-{}, {}:{}-{}",
                aln.seq1_name(),
                aln.seq1_start_from_forward(),
                aln.seq1_end_from_forward(),
                aln.seq2_name(),
                aln.seq2_start_from_forward(),
                aln.seq2_end_from_forward()
            );
        }
        let count = count_repetitiveness(&contigs, &alns);
        let (max_contig, start, end) = find_maximum(&count);
        debug!("Remove:{}:{}-{}", max_contig, start, end);
        let contig = contigs
            .iter()
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
                let seq1 = aln.seq1_name() != max_contig
                    || aln.seq1_start_from_forward() < start
                    || end < aln.seq1_end_from_forward();
                let seq2 = aln.seq2_name() != max_contig
                    || aln.seq2_start_from_forward() < start
                    || end < aln.seq2_end_from_forward();
                seq1 && seq2
            })
            .collect::<Vec<_>>();
    }
    let stdout = std::io::stdout();
    use std::io::BufWriter;
    let mut wtr = fasta::Writer::new(BufWriter::new(stdout.lock()));
    for record in contigs
        .into_iter()
        .enumerate()
        .map(|(idx, r)| {
            let seq: Vec<_> = r
                .seq()
                .iter()
                .zip(mask[r.id()].iter())
                .filter_map(|(&b, &m)| if m { Some(b) } else { None })
                .collect();
            fasta::Record::with_data(&format!("{}", idx), &None, &seq)
        })
        .chain(result.into_iter())
    {
        wtr.write_record(&record)?;
    }
    Ok(())
}

fn find_maximum<'a>(count: &'a HashMap<String, Vec<usize>>) -> (&'a str, usize, usize) {
    let (mut count_max, mut range_max, mut start, mut end, mut contig) = (0, 0, 0, 0, "");
    for (key, val) in count {
        let (cm, rm, s, e) = find_max_region(val);
        if count_max < cm || (count_max == cm && range_max < rm) {
            count_max = cm;
            range_max = rm;
            contig = key;
            start = s;
            end = e;
        }
    }
    (contig, start, end)
}
fn find_max_region(xs: &[usize]) -> (usize, usize, usize, usize) {
    let (idx, max) = xs
        .iter()
        .enumerate()
        .max_by_key(|(_, count)| *count)
        .unwrap();
    let (mut start, mut end) = (idx, idx);
    while 0 < start && 0 < xs[start] {
        start -= 1;
    }
    while end < xs.len() && 0 < xs[end] {
        end += 1;
    }
    debug!("{}", xs[idx]);
    (*max, end - start, start, end)
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
