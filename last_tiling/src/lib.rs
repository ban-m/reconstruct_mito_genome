//! This is a libray to tiling a fasta file into
//! chunks of contigs(units).
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rmp_serde;
#[macro_use]
extern crate serde;
extern crate bio_utils;
use bio_utils::fasta;
mod peak;
pub use peak::UnitDefinitions;
pub mod unit;
pub use peak::SUBUNIT_SIZE;
pub use unit::EncodedRead;
pub mod lasttab;
use lasttab::LastTAB;
use std::collections::HashMap;
use std::path::Path;
use unit::*;

/// The function to parse Last's TAB format. To see more detail,
/// see [lasttab] module.
/// # Example
/// ```ignore
/// let path = "/path/to/last/tab/file.tab";
/// let tabs = last_tiling::parse_tab_file(path).unwrap();
/// println!("{}",tabs[0].seq1_name());
/// ```
pub fn parse_tab_file<P: AsRef<Path>>(tab_file: P) -> std::io::Result<Vec<LastTAB>> {
    let lines = std::fs::read_to_string(tab_file)?;
    Ok(lines
        .lines()
        .filter(|e| !e.starts_with('#'))
        .filter_map(LastTAB::from_line)
        .collect())
}

pub fn parse_peak_file<P: AsRef<Path>>(
    peak_file: P,
    ctg_file: P,
) -> std::io::Result<UnitDefinitions> {
    let ctgs: Vec<_> = bio_utils::fasta::parse_into_vec(ctg_file)?;
    let peaks = std::fs::read_to_string(peak_file)?;
    Ok(peak::UnitDefinitions::open_peak_with_contigs(peaks, ctgs))
}

pub fn encoding(
    fasta: &[fasta::Record],
    defs: &UnitDefinitions,
    alns: &[LastTAB],
) -> Vec<EncodedRead> {
    // Distribute alignments to each reads.
    // bucket[i] is the alignment for fasta[i].
    let buckets = distribute(fasta, alns);
    debug!("There are {} buckets.", buckets.len());
    buckets
        .into_iter()
        .zip(fasta.iter())
        .map(|(bucket, seq)| into_encoding(bucket, seq, defs))
        .take(1)
        .inspect(|read| debug!("{}", read))
        .collect()
}

fn into_encoding(
    bucket: Vec<&LastTAB>,
    seq: &fasta::Record,
    defs: &UnitDefinitions,
) -> EncodedRead {
    if bucket.is_empty() {
        let read = vec![ChunkedUnit::Gap(GapUnit::new(seq.seq()))];
        return EncodedRead::from(seq.id().to_string(), read);
    }
    debug!("Encoding {} alignments", bucket.len());
    let bucket = filter_contained_alignment(bucket);
    debug!("Filter contained. Remain {} alignments", bucket.len());
    let mut start_pos = 0;
    let mut read = vec![];
    let bases = seq.seq();
    for w in bucket.windows(2) {
        // Determine whether we can useentire alignment of w[0].
        let (mut encodes, start, end) = if w[0].score() > w[1].score()
            || w[0].seq2_end_from_forward() < w[1].seq2_start_from_forward()
        {
            aln_to_encode(&w[0], w[0].seq2_end_from_forward(), defs, bases)
        } else {
            aln_to_encode(&w[1], w[1].seq2_start_from_forward(), defs, bases)
        };
        if start_pos < start {
            let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start]));
            read.push(gapunit);
        }
        read.append(&mut encodes);
        start_pos = end;
    }
    if let Some(last) = bucket.last() {
        let (mut encodes, start, end) =
            aln_to_encode(last, last.seq2_end_from_forward(), defs, bases);
        if start_pos < start {
            let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start]));
            read.push(gapunit);
        }
        read.append(&mut encodes);
        if end < bases.len() {
            let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[end..]));
            read.push(gapunit);
        }
    }
    unit::EncodedRead::from(seq.id().to_string(), read)
}

#[inline]
pub fn revcmp(seq: &[u8]) -> Vec<u8> {
    seq.into_iter()
        .rev()
        .map(|&e| match e {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

fn chop_reference_into_chunk(
    def: &UnitDefinitions,
    id: &str,
    start: usize,
    end: usize,
) -> (usize, usize, Vec<ChunkedUnit>) {
    // First, determine the location to start tiling.
    let encode_start = {
        // If the first subunit crossing boundary, skip to the next unit.
        let current_unit = def.definition(id, start).unwrap();
        let next_unit = def.definition(id, start + SUBUNIT_SIZE).unwrap();
        if current_unit == next_unit {
            current_unit.start()
                + ((start - current_unit.start()) / SUBUNIT_SIZE + 1) * SUBUNIT_SIZE
        } else {
            next_unit.start()
        }
    };
    assert!(start <= encode_start);
    let mut chunks = vec![];
    let mut pos = encode_start;
    while pos + SUBUNIT_SIZE < end {
        let unit = def.definition(id, pos).unwrap();
        assert!(pos <= unit.end(), "{:?},{}", unit, pos);
        let remaining = unit.end() - pos;
        if remaining < SUBUNIT_SIZE && end < SUBUNIT_SIZE + pos + remaining {
            // Too near the boundary and there's no additional unit.
            break;
        } else if remaining < SUBUNIT_SIZE {
            // Too near the boundary. Go to next unit.
            chunks.push(ChunkedUnit::Gap(GapUnit::new(&vec![b'-'; remaining])));
            pos = unit.end();
        } else {
            let subunit = ((pos - unit.start()) / SUBUNIT_SIZE) as u16;
            chunks.push(ChunkedUnit::En(Encode::sketch(
                unit.contig(),
                unit.num(),
                subunit,
            )));
            pos += SUBUNIT_SIZE;
        }
    }
    (encode_start, pos, chunks)
}

// Stop is the location where the tiling stops, at `seq`.
fn aln_to_encode(
    aln: &LastTAB,
    stop: usize,
    def: &UnitDefinitions,
    seq: &[u8],
) -> (Vec<ChunkedUnit>, usize, usize) {
    debug!("Refr:{}-{}", aln.seq1_start(), aln.seq1_end_from_forward());
    debug!(
        "Read:{}-{}",
        aln.seq2_start_from_forward(),
        aln.seq2_end_from_forward()
    );
    let seq = match aln.seq2_direction() {
        lasttab::Strand::Forward => seq.to_vec(),
        lasttab::Strand::Reverse => revcmp(seq),
    };
    let ctgname = aln.seq1_name();
    let refr = def.get_reference_sequence(aln.seq1_name()).unwrap().seq();
    debug!("Ctgname:{}", ctgname);
    // First, chunk the reference into subunits.
    debug!("{:?}", aln);
    let (rs, qs) = recover(aln, &refr, &seq);
    let dig = 100;
    for i in 0..rs.len() / dig {
        debug!("{}", String::from_utf8_lossy(&rs[i * dig..(i + 1) * dig]));
        debug!("{}", String::from_utf8_lossy(&qs[i * dig..(i + 1) * dig]));
        debug!("");
    }
    let (ref_encode_start, _, chunks) =
        chop_reference_into_chunk(def, ctgname, aln.seq1_start(), aln.seq1_end_from_forward());
    debug!("{:?}", chunks);
    let (mut ops, read_encode_start) = seek_to_head(aln, ref_encode_start);
    let mut read_pos = read_encode_start;
    let mut refr_pos = ref_encode_start;
    // debug!(
    //     "Read:{}",
    //     String::from_utf8_lossy(&seq[aln.seq2_start()..read_pos])
    // );
    // debug!(
    //     "Refr:{}",
    //     String::from_utf8_lossy(&refr[aln.seq1_start()..refr_pos])
    // );

    let chunks: Vec<_> = chunks
        .into_iter()
        .filter_map(|chunk| {
            debug!("Refr:{}, Read:{}", refr_pos, read_pos);
            if stop < read_pos {
                return None;
            }
            match chunk {
                ChunkedUnit::Gap(mut gu) => {
                    // The length in the reference.
                    let ref_len = gu.len();
                    //  the length in the read.
                    let (read_len, _) = seek_len(ref_len, &mut ops);
                    gu.set_bases(&seq[read_pos..read_pos + read_len]);

                    read_pos += read_len;
                    refr_pos += ref_len;
                    Some(ChunkedUnit::Gap(gu))
                }
                ChunkedUnit::En(mut encode) => {
                    let ref_len = SUBUNIT_SIZE;
                    let (read_len, operations) = seek_len(ref_len, &mut ops);
                    encode.set_bases(&seq[read_pos..read_pos + read_len]);
                    encode.set_ops(&operations);
                    read_pos += read_len;
                    refr_pos += SUBUNIT_SIZE;
                    Some(ChunkedUnit::En(encode))
                }
            }
        })
        .collect();
    (chunks, read_encode_start, read_pos)
}

// Seek ops to `len` length.
fn seek_len(len: usize, ops: &mut Vec<Op>) -> (usize, Vec<Op>) {
    debug!("{}",len);
    let mut read_len = 0;
    let mut refr_len = 0;
    let mut popped_ops = vec![];
    while refr_len < len {
        // debug!("{},{}",refr_len,read_len);
        match ops.pop().unwrap() {
            Op::Match(l) => {
                if len < refr_len + l {
                    ops.push(Op::Match(refr_len + l - len));
                    popped_ops.push(Op::Match(len - refr_len));
                    read_len += len - refr_len;
                    refr_len = len;
                } else {
                    refr_len += l;
                    read_len += l;
                    popped_ops.push(Op::Match(l));
                }
            }
            Op::Seq1In(l) => {
                read_len += l;
                popped_ops.push(Op::Seq1In(l));
            }
            Op::Seq2In(l) => {
                if len < refr_len + l {
                    ops.push(Op::Seq2In(refr_len + l - len));
                    popped_ops.push(Op::Seq2In(len - refr_len));
                    refr_len = len;
                } else {
                    refr_len += l;
                    popped_ops.push(Op::Seq2In(l));
                }
            }
        }
    }
    (read_len, popped_ops)
}

use lasttab::Op;
fn seek_to_head(aln: &LastTAB, ref_encode_start: usize) -> (Vec<Op>, usize) {
    let mut ops = aln.alignment();
    ops.reverse();
    let refr_pos = aln.seq1_start();
    let (len, _) = seek_len(ref_encode_start - refr_pos, &mut ops);
    (ops, aln.seq2_start() + len)
}

fn filter_contained_alignment<'a>(mut bucket: Vec<&'a LastTAB>) -> Vec<&'a LastTAB> {
    bucket.sort_by_key(|aln| aln.seq2_start_from_forward());
    let (mut start, mut end) = (0, 0);
    bucket
        .into_iter()
        .filter(|aln| {
            let s = aln.seq2_start_from_forward();
            let e = aln.seq2_start_from_forward();
            assert!(start <= s);
            if e < end {
                start = s;
                false
            } else {
                start = s;
                end = e;
                true
            }
        })
        .collect()
}

fn distribute<'a>(fasta: &[fasta::Record], alns: &'a [LastTAB]) -> Vec<Vec<&'a LastTAB>> {
    let mut alignments_bucket: Vec<Vec<&LastTAB>> = vec![vec![]; fasta.len()];
    let id_to_idx: HashMap<_, _> = fasta
        .iter()
        .map(|e| e.id())
        .enumerate()
        .map(|(idx, id)| (id, idx))
        .collect();
    for aln in alns
        .iter()
        .filter(|aln| aln.alignment_length() > SUBUNIT_SIZE)
    {
        alignments_bucket[id_to_idx[aln.seq2_name()]].push(aln);
    }
    alignments_bucket
}
fn recover(aln: &LastTAB, refr: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let (mut r, mut q) = (aln.seq1_start(), aln.seq2_start());
    let ops = aln.alignment();
    let (mut rs, mut qs) = (vec![], vec![]);
    for op in ops {
        match op {
            Op::Match(l) => {
                rs.extend(&refr[r..r + l]);
                qs.extend(&query[q..q + l]);
                r += l;
                q += l;
            }
            Op::Seq1In(l) => {
                rs.extend(&vec![b'-'; l]);
                qs.extend(&query[q..q + l]);
                q += l;
            }
            Op::Seq2In(l) => {
                qs.extend(&vec![b'-'; l]);
                rs.extend(&refr[r..r + l]);
                r += l;
            }
        }
    }
    (rs, qs)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
    #[test]
    fn seek_test() {
        use lasttab::Op::*;
        let ops = vec![
            Match(10),
            Seq1In(3),
            Seq2In(4),
            Match(10),
            Seq1In(5),
            Match(2),
            Seq2In(3),
            Match(4),
        ];
        let len = 4;
        let mut res = ops.clone();
        let (r_len, popped) = seek_len(len, &mut res);
        assert_eq!(r_len, 4);
        assert_eq!(popped, vec![Match(4)]);
        {
            let mut ans = ops.clone();
            ans.pop();
            assert_eq!(res, ans);
        }
        let len = 14;
        let mut res = ops.clone();
        let (r_len, popped) = seek_len(len, &mut res);
        assert_eq!(r_len, 11);
        assert_eq!(
            popped,
            vec![Match(4), Seq2In(3), Match(2), Seq1In(5), Match(5)]
        );
        {
            let mut ans = ops.clone();
            (0..5).for_each(|_| {
                ans.pop().unwrap();
            });
            ans.push(Match(5));
            assert_eq!(res, ans);
        }
        let len = 5;
        let mut res = ops.clone();
        let (r_len, popped) = seek_len(len, &mut res);
        assert_eq!(r_len, 4);
        assert_eq!(popped, vec![Match(4), Seq2In(1)]);
        {
            let mut ans = ops.clone();
            ans.pop();
            ans.pop();
            ans.push(Seq2In(2));
            assert_eq!(res, ans);
        }
    }
}
