//! This is a libray to tiling a fasta file into
//! chunks of contigs(units).
#[allow(unused_imports)]
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
pub const UNIT_SIZE: usize = 150;
// If an alignment is in a repetitive region,
// we check whether the read also aligns at CHECK_POINT up stream.
// If it does not, we discard that alignment. Otherwise, we keep it.
const CHECK_POINT: usize = 550;
/// Parse given aln files into repeats. It needs contig information such as index of a contig.
const THR: usize = 1_000;
pub fn into_repeats(alns: &[LastTAB], contig: &Contigs) -> Vec<RepeatPairs> {
    alns.into_iter()
        .filter(|aln| {
            // Check complete alignment.
            let seq1_cmp = aln.seq1_matchlen() != aln.seq1_len();
            let seq2_cmp = aln.seq2_matchlen() != aln.seq2_len();
            // Check Long alignment.
            let long = aln.seq1_matchlen() > THR && aln.seq2_matchlen() > THR;
            seq1_cmp && seq2_cmp && long
        })
        .filter_map(|aln| repeat::RepeatPairs::new(&aln, &contig))
        .collect()
}

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

pub fn encoding(fasta: &[fasta::Record], defs: &Contigs, alns: &[LastTAB]) -> Vec<EncodedRead> {
    // Distribute alignments to each reads.
    // bucket[i] is the alignment for fasta[i].
    let buckets = distribute(fasta, alns);
    // debug!("There are {} buckets.", buckets.len());
    let buckets: Vec<_> = buckets.into_iter().zip(fasta.iter()).collect();
    buckets
        .into_iter()
        .map(|(bucket, seq)| {
            if bucket.is_empty() {
                let read = vec![ChunkedUnit::Gap(GapUnit::new(seq.seq(), None))];
                EncodedRead::from(seq.id().to_string(), read)
            } else {
                into_encoding(bucket, seq, defs)
            }
        })
        .collect()
}

pub fn encoding_w_repeat(
    fasta: &[fasta::Record],
    defs: &Contigs,
    alns: &[LastTAB],
    repeat: &[RepeatPairs],
) -> Vec<EncodedRead> {
    // Distribute alignments to each reads.
    // bucket[i] is the alignment for fasta[i].
    let buckets = distribute(fasta, alns);
    // debug!("There are {} buckets.", buckets.len());
    let buckets: Vec<_> = buckets.into_iter().zip(fasta.iter()).collect();
    buckets
        .into_iter()
        .map(|(bucket, seq)| {
            let bucket = trim_aln_in_repetitive(bucket, repeat);
            if bucket.is_empty() {
                let read = vec![ChunkedUnit::Gap(GapUnit::new(seq.seq(), None))];
                EncodedRead::from(seq.id().to_string(), read)
            } else {
                into_encoding(bucket, seq, defs)
            }
        })
        .collect()
}

fn trim_aln_in_repetitive<'a>(
    bucket: Vec<&'a LastTAB>,
    repeat: &[RepeatPairs],
) -> Vec<&'a LastTAB> {
    bucket
        .iter()
        .filter_map(|aln| match is_in_repeat(aln, repeat) {
            None => Some(aln),
            Some(rep) if has_flanking(rep, &bucket) => Some(aln),
            Some(_) => None,
        })
        .copied()
        .collect()
}
use repeat::Repeat;
fn is_in_repeat<'a>(aln: &LastTAB, repeat: &'a [RepeatPairs]) -> Option<&'a Repeat> {
    let margin = 200;
    let (start, end) = (aln.seq1_start_from_forward(), aln.seq1_end_from_forward());
    let contig = aln.seq1_name();
    //debug!("Checking {}-{}", start, end);
    repeat
        .iter()
        .filter_map(|rs| {
            rs.inner()
                .iter()
                .filter(|r| {
                    let r_start = r.start().max(margin) - margin;
                    let r_end = r.end() + margin;
                    let b = r.name() == contig && r_start <= start && end <= r_end;
                    b
                })
                .nth(0)
        })
        .nth(0)
}

fn has_flanking(rep: &Repeat, bucket: &[&LastTAB]) -> bool {
    let check_start = rep.start().max(CHECK_POINT) - CHECK_POINT;
    let check2_start = check_start.max(CHECK_POINT) - CHECK_POINT;
    let count = bucket
        .iter()
        .filter(|aln| {
            let (start, end) = (aln.seq1_start_from_forward(), aln.seq1_end_from_forward());
            let c1 = start < check_start && check_start < end;
            let c2 = start < check2_start && check2_start < end;
            c1 || c2
        })
        .count();
    count != 0
}

fn into_encoding(bucket: Vec<&LastTAB>, seq: &fasta::Record, defs: &Contigs) -> EncodedRead {
    // debug!("Encoding {} alignments", bucket.len());
    // debug!("Read:{},{}len", seq.id(), seq.seq().len());
    let bucket = filter_contained_alignment(bucket, defs);
    // debug!("Filter contained. Remain {} alignments", bucket.len());
    // for aln in &bucket {
    //     debug!(
    //         "{}-{}({}:{}-{})",
    //         aln.seq2_start_from_forward(),
    //         aln.seq2_end_from_forward(),
    //         aln.seq1_name(),
    //         aln.seq1_start_from_forward(),
    //         aln.seq1_end_from_forward()
    //     );
    // }
    let (mut start_pos, mut read) = (0, vec![]);
    let bases = seq.seq();
    // Id of bucket[buclet.len()-2]. If bucket.len()==1, it should be None.
    let mut prev_contig_id = defs.get_id(bucket[0].seq1_name()).unwrap();
    // The procedure here would be a little bit tricky.
    let mut bucket = bucket.iter();
    let mut target = bucket.next().unwrap();
    while let Some(next) = bucket.next() {
        let (sp, units) =
            encoding_single_alignment(target, next, defs, bases, start_pos, prev_contig_id);
        read.extend(units);
        start_pos = sp;
        target = next;
        prev_contig_id = defs.get_id(target.seq1_name()).unwrap();
    }
    let last = target;
    let (encodes, start, end) = aln_to_encode(last, last.seq2_end_from_forward(), defs, bases);
    let c = defs.get_id(&last.seq1_name()).unwrap();
    // debug!("Start from {} in read", start);
    let (start, encodes) = pop_encodes_until(start, start_pos, encodes);
    if start_pos < start {
        let cs = (c.min(prev_contig_id), c.max(prev_contig_id));
        let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start], Some(cs)));
        read.push(gapunit);
    }
    read.extend(encodes);
    // debug!("SP:{}->{}", start_pos, end.max(start_pos));
    start_pos = end.max(start_pos);
    if start_pos < bases.len() {
        let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..], Some((c, c))));
        read.push(gapunit);
    }
    unit::EncodedRead::from(seq.id().to_string(), read)
}

type ChunkedUnits = Vec<ChunkedUnit>;
fn encoding_single_alignment(
    target: &LastTAB,
    next: &LastTAB,
    defs: &Contigs,
    bases: &[u8],
    start_pos: usize,
    prev: u16,
) -> (usize, ChunkedUnits) {
    assert!(!check_contained(target, next));
    let former_stop = target.seq2_end_from_forward();
    let later_start = next.seq2_start_from_forward();
    let (encodes, start, end) = if former_stop < later_start {
        // No overlap.
        aln_to_encode(target, target.seq2_end_from_forward(), defs, bases)
    } else if defs.get_id(target.seq1_name()).unwrap() < defs.get_id(next.seq1_name()).unwrap() {
        // Overlap.Take the younger contig.
        aln_to_encode(target, target.seq2_end_from_forward(), defs, bases)
    } else {
        aln_to_encode(target, next.seq2_start_from_forward(), defs, bases)
    };
    // If there are overlapping, there would be some "doubly encoded" regions.
    let (start, mut encodes) = pop_encodes_until(start, start_pos, encodes);
    // debug!("Start from {} in read", start);
    let mut units = vec![];
    if start_pos < start {
        let cs = defs
            .get_id(target.seq1_name())
            .map(|c| (c.min(prev), c.max(prev)));
        let gapunit = ChunkedUnit::Gap(GapUnit::new(&bases[start_pos..start], cs));
        units.push(gapunit);
    }
    units.append(&mut encodes);
    // debug!("SP:{}->{}", start_pos, end.max(start_pos));
    (end.max(start_pos), units)
}

fn pop_encodes_until(
    mut encode_start: usize,
    actual_start: usize,
    encodes: ChunkedUnits,
) -> (usize, ChunkedUnits) {
    let mut encodes = encodes.into_iter();
    while encode_start < actual_start {
        match encodes.next() {
            Some(unit) => encode_start += unit.len(),
            None => break,
        }
    }
    (encode_start, encodes.collect())
}

#[allow(dead_code)]
fn check_contained(a1: &LastTAB, a2: &LastTAB) -> bool {
    // Check whether a1 covers a2.
    let (s, t) = (a1.seq2_start_from_forward(), a1.seq2_end_from_forward());
    let (x, y) = (a2.seq2_start_from_forward(), a2.seq2_end_from_forward());
    (s < x) && (y < t)
}
#[allow(dead_code)]
fn is_not_contained_by(rs: &RepeatPairs, a: &LastTAB, s: usize, t: usize) -> bool {
    rs.inner().iter().all(|r| is_not_contained(r, a, s, t))
}

#[allow(dead_code)]
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
    // debug!("Refr:{}-{}", aln.seq1_start(), aln.seq1_end_from_forward());
    // debug!(
    //     "Read:{}-{}",
    //     aln.seq2_start_from_forward(),
    //     aln.seq2_end_from_forward()
    // );
    // let ctgname = aln.seq1_name();
    // Flip reference rather than query.
    // let _refr = match aln.seq2_direction() {
    //     lasttab::Strand::Forward => def.get(ctgname).unwrap(),
    //     lasttab::Strand::Reverse => def.get_revcmp(ctgname).unwrap(),
    // };
    // debug!("Ctgname:{}", ctgname);
    // First, chunk the reference into subunits.
    let (ref_encode_start, _, chunks) = chop_reference_into_chunk(def, aln);
    // debug!("{:?},{}", chunks, ref_encode_start);
    if chunks.is_empty() {
        return (vec![], 0, stop);
    }
    let (mut ops, read_encode_start) = seek_to_head(aln, ref_encode_start);
    let mut read_pos = read_encode_start;
    let mut refr_pos = ref_encode_start;
    {
        // let start = ref_start_end(aln).0;
        // debug!(
        //     "Read:{}",
        //     String::from_utf8_lossy(&seq[aln.seq2_start_from_forward()..read_pos])
        // );
        // debug!("Refr:{}", String::from_utf8_lossy(&refr[start..refr_pos]));
    }
    let chunks: Vec<_> = chunks.into_iter().fold(vec![], |mut chunks, mut encode| {
        // debug!("Refr:{}, Read:{}", refr_pos, read_pos);
        let ref_len = UNIT_SIZE;
        let (read_len, operations) = seek_len(ref_len, &mut ops);
        if read_pos + read_len < stop {
            encode.set_bases(&seq[read_pos..read_pos + read_len]);
            encode.set_ops(&operations); // This is the vectrized version of operations. Not stack-version.
            refr_pos += UNIT_SIZE;
            read_pos += read_len;
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

fn filter_contained_alignment<'a>(
    mut bucket: Vec<&'a LastTAB>,
    _defs: &Contigs,
) -> Vec<&'a LastTAB> {
    use std::cmp::Ordering;
    bucket.sort_by(|aln1, aln2| {
        if aln1.seq2_start_from_forward() < aln2.seq2_start_from_forward() {
            Ordering::Less
        } else if aln1.seq2_start_from_forward() > aln2.seq2_start_from_forward() {
            Ordering::Greater
        } else if aln1.seq2_end_from_forward() > aln2.seq2_end_from_forward() {
            Ordering::Less
        } else if aln1.seq2_end_from_forward() < aln2.seq2_end_from_forward() {
            Ordering::Greater
        } else if aln1.score() > aln2.score() {
            Ordering::Less
        } else if aln1.score() < aln2.score() {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });
    let mut end = 0;
    bucket
        .into_iter()
        .filter(|aln| {
            let e = aln.seq2_end_from_forward();
            // let order = defs.get_id(aln.seq1_name()).unwrap();
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
        if let Some(idx) = id_to_idx.get(aln.seq2_name()) {
            alignments_bucket[*idx].push(aln);
        }
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
        assert_eq!(r_len, 16);
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
