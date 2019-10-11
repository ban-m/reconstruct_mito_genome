//! This is a libray to tiling a fasta file into
//! chunks of contigs(units).
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rmp_serde;
#[macro_use]
extern crate serde;
extern crate bio_utils;
extern crate rayon;
pub mod contig;
pub mod lasttab;
pub mod repeat;
pub mod unit;

pub use contig::Contigs;
pub use lasttab::LastTAB;
pub use lasttab::Op;
pub use unit::EncodedRead;

use bio_utils::fasta;
use rayon::prelude::*;
use repeat::RepeatPairs;
use std::collections::HashMap;
use std::path::Path;
use unit::*;
pub const UNIT_SIZE: usize = 200;

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

pub fn remove_repeats(alns: Vec<LastTAB>, defs: &Contigs, rep: &[RepeatPairs]) -> Vec<LastTAB> {
    let repeats: HashMap<u16, Vec<(usize, usize)>> =
        rep.iter().fold(HashMap::new(), |mut map, reps| {
            reps.inner().iter().for_each(|rep| {
                let entry = map.entry(rep.id()).or_default();
                entry.push((rep.start(), rep.end()));
            });
            map
        });
    alns.into_par_iter()
        .filter(|aln| {
            let id = match defs.get_id(aln.seq1_name()) {
                Some(res) => res,
                None => {
                    eprintln!("Error. Seq1 of {} is not in the referneces.", aln);
                    return false;
                }
            };
            let start = aln.seq1_start_from_forward();
            let end = aln.seq1_end_from_forward();
            match repeats.get(&id) {
                Some(ref reps) => !reps.iter().any(|&(s, e)| s < start && end < e),
                // Corresponds to "no repeats there."
                None => true,
            }
        })
        .collect()
}

pub fn encoding(
    fasta: &[fasta::Record],
    defs: &Contigs,
    alns: &[LastTAB],
    repeats: &[RepeatPairs],
) -> Vec<EncodedRead> {
    // Distribute alignments to each reads.
    // bucket[i] is the alignment for fasta[i].
    let buckets = distribute(fasta, alns);
    debug!("There are {} buckets.", buckets.len());
    let buckets: Vec<_> = buckets.into_iter().zip(fasta.iter()).collect();
    buckets
        .into_par_iter()
        .map(|(bucket, seq)| into_encoding(bucket, seq, defs, repeats))
        // .enumerate()
        // .inspect(|(idx, read)| debug!("{},{}", idx, read))
        // .map(|(_, read)| read)
        .collect()
}

fn into_encoding(
    bucket: Vec<&LastTAB>,
    seq: &fasta::Record,
    defs: &Contigs,
    repeats: &[RepeatPairs],
) -> EncodedRead {
    if bucket.is_empty() {
        let read = vec![ChunkedUnit::Gap(GapUnit::new(seq.seq()))];
        return EncodedRead::from(seq.id().to_string(), read);
    }
    debug!("Encoding {} alignments", bucket.len());
    debug!("Read:{},{}len", seq.id(), seq.seq().len());
    let bucket = filter_contained_alignment(bucket);
    debug!("Filter contained. Remain {} alignments", bucket.len());
    for aln in &bucket {
        debug!(
            "{}-{}({}:{}-{})",
            aln.seq2_start_from_forward(),
            aln.seq2_end_from_forward(),
            aln.seq1_name(),
            aln.seq1_start_from_forward(),
            aln.seq1_end_from_forward()
        );
    }
    let (mut start_pos, mut read) = (0, vec![]);
    let bases = seq.seq();
    for w in bucket.windows(2) {
        // Determine whether we can use entire alignment of w[0].
        let (former_start, former_stop) = ref_start_end(&w[0]);
        let (later_start, later_stop) = ref_start_end(&w[1]);
        let (mut encodes, start, end) = if former_stop < later_start {
            // No overlap.
            aln_to_encode(&w[0], w[0].seq2_end_from_forward(), defs, bases)
        } else {
            if repeats
                .iter()
                .all(|rp| is_not_contained_by(rp, &w[0], former_start, former_stop))
            {
                aln_to_encode(&w[0], w[0].seq2_end_from_forward(), defs, bases)
            } else {
                aln_to_encode(&w[1], w[1].seq2_start_from_forward(), defs, bases)
            }
        };
        if start_pos < start {
            let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start]));
            read.push(gapunit);
        }
        read.append(&mut encodes);
        debug!("SP:{}->{}", start_pos, end);
        start_pos = end;
    }
    let last = bucket.last().unwrap();
    let (mut encodes, start, end) = aln_to_encode(last, last.seq2_end_from_forward(), defs, bases);
    if start_pos < start {
        let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start]));
        read.push(gapunit);
    }
    read.append(&mut encodes);
    if end < bases.len() {
        let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[end..]));
        read.push(gapunit);
    }
    unit::EncodedRead::from(seq.id().to_string(), read)
}

fn is_not_contained_by(rs: &RepeatPairs, a: &LastTAB, s: usize, t: usize) -> bool {
    rs.inner().iter().all(|r| is_not_contained(r, a, s, t))
}

#[inline]
fn is_not_contained(r: &repeat::Repeat, a: &LastTAB, s: usize, t: usize) -> bool {
    r.name() != a.seq2_name() || s < r.start() || r.end() < t
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

// (ChunkedAlignment, Start position of encoding, End position of encoding).
// These locations are the position at the read, not references.
type AlnToEncode = (Vec<ChunkedUnit>, usize, usize);
// Stop is the location where the tiling stops, at `seq`.
fn aln_to_encode(aln: &LastTAB, stop: usize, def: &Contigs, seq: &[u8]) -> AlnToEncode {
    debug!("Refr:{}-{}", aln.seq1_start(), aln.seq1_end_from_forward());
    debug!(
        "Read:{}-{}",
        aln.seq2_start_from_forward(),
        aln.seq2_end_from_forward()
    );
    let ctgname = aln.seq1_name();
    // Flip reference rather than query.
    let refr = match aln.seq2_direction() {
        lasttab::Strand::Forward => def.get(ctgname).unwrap(),
        lasttab::Strand::Reverse => def.get_revcmp(ctgname).unwrap(),
    };
    debug!("Ctgname:{}", ctgname);
    // First, chunk the reference into subunits.
    let (ref_encode_start, _, chunks) = chop_reference_into_chunk(def, aln);
    debug!("{:?},{}", chunks, ref_encode_start);
    if chunks.is_empty() {
        return (vec![], 0, stop);
    }
    let (mut ops, read_encode_start) = seek_to_head(aln, ref_encode_start);
    let mut read_pos = read_encode_start;
    let mut refr_pos = ref_encode_start;
    {
        let start = ref_start_end(aln).0;
        debug!(
            "Read:{}",
            String::from_utf8_lossy(&seq[aln.seq2_start_from_forward()..read_pos])
        );
        debug!("Refr:{}", String::from_utf8_lossy(&refr[start..refr_pos]));
    }
    let chunks: Vec<_> = chunks.into_iter().fold(vec![], |mut chunks, mut encode| {
        // debug!("Refr:{}, Read:{}", refr_pos, read_pos);
        if read_pos < stop {
            let ref_len = UNIT_SIZE;
            let (read_len, operations) = seek_len(ref_len, &mut ops);
            encode.set_bases(&seq[read_pos..read_pos + read_len]);
            encode.set_ops(&operations); // This is the vectrized version of operations. Not stack-version.
            read_pos += read_len;
            refr_pos += UNIT_SIZE;
            chunks.push(ChunkedUnit::En(encode));
        }
        chunks
    });
    (chunks, read_encode_start, read_pos)
}

// Convert alignment into array of unit.
// Reference's start and end position also returned.
fn chop_reference_into_chunk(def: &Contigs, aln: &LastTAB) -> (usize, usize, Vec<Encode>) {
    // First, determine the location to start tiling.
    let (start, end) = (aln.seq1_start(), aln.seq1_end_from_forward());
    let name = aln.seq1_name();
    let id = def.get_id(name).unwrap();
    let encode_start = if start == 0 { 0 } else { start / UNIT_SIZE + 1 };
    let encode_end = end / UNIT_SIZE;
    assert!(start <= UNIT_SIZE * encode_start);
    if aln.seq2_direction().is_forward() {
        let chunks: Vec<_> = (encode_start..encode_end)
            .map(|i| Encode::sketch(id, i as u16, true))
            .collect();
        (encode_start * UNIT_SIZE, encode_end * UNIT_SIZE, chunks)
    } else {
        let len = aln.seq1_len();
        let chunks: Vec<_> = (encode_start..encode_end)
            .map(|i| Encode::sketch(id, i as u16, false))
            .rev()
            .collect();
        let start = len - encode_end * UNIT_SIZE;
        let end = len - encode_start * UNIT_SIZE;
        (start, end, chunks)
    }
}

// Seek the alignment operation up to ref_encode_start position.
// Note that if the alignment is reverse complement, the operations should be treated as such.
// The second returned value is the position of the read, which the first op would be applied.
fn seek_to_head(aln: &LastTAB, ref_encode_start: usize) -> (Vec<Op>, usize) {
    let mut ops = match aln.seq2_direction() {
        lasttab::Strand::Forward => aln.alignment(),
        lasttab::Strand::Reverse => {
            let mut a = aln.alignment();
            a.reverse();
            a
        }
    };
    // Additional reverse, which should be needed because we want to
    // treat the alignment operations like a stack, rather than a vector.
    ops.reverse();
    let start = ref_start_end(aln).0;
    let (len, _) = seek_len(ref_encode_start - start, &mut ops);
    (ops, aln.seq2_start_from_forward() + len)
}

// [start..end) of reference. If the alignment is reverse complement,
// the coordinate would be reversed.
fn ref_start_end(aln: &LastTAB) -> (usize, usize) {
    match aln.seq2_direction() {
        lasttab::Strand::Forward => (aln.seq1_start(), aln.seq1_end_from_forward()),
        lasttab::Strand::Reverse => {
            let len = aln.seq1_len();
            let mlen = aln.seq1_matchlen();
            let start = len - mlen - aln.seq1_start();
            (start, start + mlen)
        }
    }
}

// Seek ops to `len` length in the reference. Here, one need not to consider
// whether the alignment is reverse or forward. It should be
// treated beforehand(such as in `seek_head`).
// The first returned value is the consumed length of the query.
// The second is the consumed operations.
// Input: [Op4 Op3 Op2 Op1] <- From this direction the operation should be applied.
// Output: [Op1 Op2 Op3'] <- The last popped operation is the last element of the retuned value.
// The remaining ops is [Op4 Op3']. Here, Op3 is splited into two operations(every ops has its length, men).
fn seek_len(len: usize, ops: &mut Vec<Op>) -> (usize, Vec<Op>) {
    // debug!("{}",len);
    let mut read_len = 0;
    let mut refr_len = 0;
    let mut popped_ops = vec![];
    while refr_len < len {
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

fn filter_contained_alignment<'a>(mut bucket: Vec<&'a LastTAB>) -> Vec<&'a LastTAB> {
    use std::cmp::Ordering;
    bucket.sort_by(|aln1, aln2| {
        if aln1.seq2_start_from_forward() < aln2.seq2_start_from_forward() {
            Ordering::Less
        } else if aln1.seq2_start_from_forward() > aln2.seq2_start_from_forward() {
            Ordering::Greater
        } else if aln1.score() > aln2.score() {
            Ordering::Less
        } else if aln1.score() < aln2.score() {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });
    for aln in &bucket {
        debug!(
            "{}-{}({}:{}-{})",
            aln.seq2_start_from_forward(),
            aln.seq2_end_from_forward(),
            aln.seq1_name(),
            aln.seq1_start_from_forward(),
            aln.seq1_end_from_forward()
        );
    }
    let mut end = 0;
    bucket
        .into_iter()
        .filter(|aln| {
            let e = aln.seq2_end_from_forward();
            if e <= end + 1 {
                false
            } else {
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
        .filter(|aln| aln.alignment_length() > 3 * UNIT_SIZE)
    {
        alignments_bucket[id_to_idx[aln.seq2_name()]].push(aln);
    }
    alignments_bucket
}

pub fn recover(aln: &LastTAB, refr: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let (mut r, mut q) = (aln.seq1_start(), aln.seq2_start());
    let ops = aln.alignment();
    let (mut rs, mut os, mut qs) = (vec![], vec![], vec![]);
    for op in ops {
        match op {
            Op::Match(l) => {
                rs.extend(&refr[r..r + l]);
                qs.extend(&query[q..q + l]);
                let o = refr[r..r + l]
                    .iter()
                    .zip(query[q..q + l].iter())
                    .map(|(a, b)| if a == b { b'|' } else { b'X' });
                os.extend(o);
                r += l;
                q += l;
            }
            Op::Seq1In(l) => {
                rs.extend(&vec![b'-'; l]);
                qs.extend(&query[q..q + l]);
                os.extend(&vec![b' '; l]);
                q += l;
            }
            Op::Seq2In(l) => {
                qs.extend(&vec![b'-'; l]);
                rs.extend(&refr[r..r + l]);
                os.extend(&vec![b' '; l]);
                r += l;
            }
        }
    }
    (rs, os, qs)
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
