extern crate bio;
extern crate env_logger;
extern crate metrohash;
extern crate run_length_encoding;
#[macro_use]
extern crate log;
use metrohash::MetroHashMap;
use run_length_encoding::BaseWithLength;
use run_length_encoding::U8Base;

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Kmer {
    pub seq: Vec<u8>,
}

impl std::fmt::Debug for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for b in self.seq.iter() {
            write!(f, "{}", *b as char)?;
        }
        Ok(())
    }
}

impl std::fmt::Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut buf = String::new();
        for b in self.seq.iter() {
            buf.push(*b as char);
        }
        write!(f, "{}", buf)
    }
}

// impl std::hash::Hash for Kmer {
//     fn hash<H>(&self, state: &mut H)
//     where
//         H: std::hash::Hasher,
//     {
//         let seq: Vec<u8> = self.seq.iter().map(U8Base::base).collect();
//         state.write(&seq);
//     }
// }

// impl std::cmp::PartialEq for Kmer {
//     fn eq(&self, other: &Self) -> bool {
//         self.len() == other.len()
//             && self
//                 .seq
//                 .iter()
//                 .zip(other.seq.iter())
//                 .all(|(x, y)| x.base() == y.base())
//     }
// }

// impl std::cmp::Eq for Kmer {}

fn cmpl(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => 0,
    }
}

// Determin whether or not the `seq` is canonical,
// in other words, return `seq < revcmp(seq)`
fn is_canonical(seq: &[u8]) -> bool {
    let len = seq.len() - 1;
    let mut idx = 0;
    loop {
        match seq[idx].cmp(&cmpl(seq[len - idx])) {
            std::cmp::Ordering::Less => break true,
            std::cmp::Ordering::Greater => break false,
            std::cmp::Ordering::Equal => idx += 1,
        }
        if idx > len / 2 {
            break false;
        };
    }
}

// Canonicalize the given kmer.
fn canonicalize(kmer: &[u8]) -> Vec<u8> {
    if is_canonical(kmer) {
        kmer.to_vec()
    } else {
        kmer.iter().rev().map(|&e| cmpl(e)).collect()
    }
}

impl Kmer {
    pub fn new(kmer: &[U8Base]) -> Self {
        Self {
            seq: canonicalize(kmer),
        }
    }
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

fn encoding(records: &[&[u8]]) -> Vec<run_length_encoding::DefaultSeq> {
    records
        .iter()
        .enumerate()
        .map(|(idx, seq)| run_length_encoding::DefaultSeq::new(idx, seq))
        .collect()
}

// Records are reference to vector of reference to sequence.
// Fisrt, it convert given records to run-length-encoded, then,
// It create hashmap with canonical k-mer.
pub fn kmer_counting(records: &[&[u8]], k: usize) -> MetroHashMap<Kmer, u16> {
    let runlength_encoded = encoding(records);
    let reads: Vec<Vec<u8>> = runlength_encoded
        .into_iter()
        .map(|e| e.seq().iter().map(U8Base::base).collect())
        .collect();
    debug!(
        "Encoded reads:{} bases",
        reads.iter().map(|e| e.len()).sum::<usize>()
    );
    let hashmap = reads
        .into_iter()
        .fold(MetroHashMap::default(), |mut res, seq| {
            seq.windows(k).map(Kmer::new).for_each(|kmer| {
                let x = res.entry(kmer).or_default();
                *x += 1;
            });
            res
        });
    debug!("Constructed Hashmap");
    hashmap
}

pub fn kmer_counting_native(records: &[Vec<u8>], k: usize) -> MetroHashMap<Kmer, u16> {
    let hashmap = records
        .into_iter()
        .fold(MetroHashMap::default(), |mut res, seq| {
            seq.windows(k).map(Kmer::new).for_each(|kmer| {
                let x = res.entry(kmer).or_default();
                *x += 1;
            });
            res
        });
    debug!("Constructed Hashmap");
    hashmap
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
