//! This is a implementation of de Bruijn Graph. It uses canonical k-mers.
//! In other words, a k-mer and its reverse complement are regarded as the same.
//!
//! There are two types of dBG, one is for exact dBG, and one is for a guided dBG.
//! More precisely, exact a dBG supports following operations:
//!
//! - Membership: k-mer -> u8 : check whether or not a dBG contains a given k-mer.
//! It supports weight. In other words, it returns the occurrence so far.
//! - In-Neighbor: k-mer -> Vec<u8>: if a char c is in the returned value, then c + x[..k-1] is
//! in the dataset.
//! - Out-Neighbor: k-mer -> Vec<u8>: Same as above.
//! - Enumerate: enumerating all the k-mers inside the structure.
//!
//! a guided dBG, on the other hand, supports following operations:
//! - In-Neighbor: k-mer -> Vec<u8>: if a char c is in the returned value, then c + x[..k-1] is
//! in the dataset.
//! - Out-Neighbor: k-mer -> Vec<u8>: Same as above.
//!
//! Both exact and guided dBG support the following operations:
//!
//! - Setup: (Int,Alphabet) -> () : set up this dBG.
//! - Insert: k-mer -> (): inserting a k-mer.
//!
extern crate metrohash;
use metrohash::MetroHashMap;

/// The trait to represent exact de Bruijn graph.
pub trait ExactdeBruijnGraph {
    /// Check the existance. Return the number of occurence.
    fn num(&self, kmer: &[u8]) -> u8;
    /// Get the in-neighbor(in other words, the kmer having the first k-1 mer of the queried k-mer as suffix.
    fn get_in_exact(&self, kmer: &[u8]) -> Vec<u8>;
    /// Get the in-neighbor(in other words, the kmer having the last k-1 mer of the queried k-mer as prefix.
    fn get_out_exact(&self, kmer: &[u8]) -> Vec<u8>;
    /// Enumerate all the k-mers in the dataset.
    fn new(k: usize, alphabet: &[u8]) -> Self;
    /// Inserting the element. It returns the occurence of the inserted k-mer so far.
    fn insert(&mut self, kmer: &[u8]) -> Option<u8>;
}

pub trait GuideddeBuruijnGraph {
    fn guide_in(&self, kmer: &[u8]) -> Vec<u8>;
    fn guide_out(&self, kmer: &[u8]) -> Vec<u8>;
    fn new(k: usize, alphabet: &[u8]) -> Self;
    fn insert(&mut self, kmer: &[u8]) -> Option<()>;
}

fn cmpl(base:u8) -> u8 {
    match &base {
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
    if seq.len() % 2 == 0 {
        eprintln!("Do not use canonicality when k is even.");
    }
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
        kmer.iter().rev().map(|&e|cmpl(e)).collect()
    }
}
#[derive(Debug)]
pub struct NaivedeBruijnGraph {
    k: usize,
    alphabet: Vec<u8>,
    map: MetroHashMap<Vec<u8>, u8>,
}

impl NaivedeBruijnGraph {
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    pub fn iter(&self)-> std::collections::hash_map::Iter<Vec<u8>,u8>{
        self.map.iter()
    }
    pub fn keys(&self) -> std::collections::hash_map::Keys<Vec<u8>,u8>{
        self.map.keys()
    }
}

impl ExactdeBruijnGraph for NaivedeBruijnGraph {
    fn new(k: usize, alphabet: &[u8]) -> Self {
        let map = MetroHashMap::default();
        Self {
            k,
            alphabet: alphabet.to_vec(),
            map,
        }
    }
    fn insert(&mut self, kmer: &[u8]) -> Option<u8> {
        if kmer.len() == self.k {
            let x = self.map.entry(canonicalize(kmer)).or_insert(0);
            *x += 1;
            Some(*x)
        } else {
            None
        }
    }
    fn num(&self, kmer: &[u8]) -> u8 {
        if is_canonical(kmer) {
            match self.map.get(kmer) {
                Some(res) => *res,
                None => 0,
            }
        } else {
            match self.map.get(&canonicalize(kmer)) {
                Some(res) => *res,
                None => 0,
            }
        }
    }
    fn get_in_exact(&self, kmer: &[u8]) -> Vec<u8> {
        self.alphabet
            .iter()
            .copied()
            .filter(|a| self.num(&vec![vec![*a], kmer[1..].to_vec()].concat()) != 0)
            .collect()
    }
    fn get_out_exact(&self, kmer: &[u8]) -> Vec<u8> {
        self.alphabet
            .iter()
            .copied()
            .filter(|a| self.num(&vec![kmer[..self.k - 1].to_vec(), vec![*a]].concat()) != 0)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
    #[test]
    fn check_naive_dbg() {
        let k = 11;
        let kmers: Vec<_> = vec![
            b"AGTACGTACGT",
            b"CAGTGTGCATG",
            b"TACGCTTCTGG",
            b"CGCCGCGCTTA",
            b"AAAAAAAAAAA",
            b"CCCCCGGGGGG",
        ]
        .into_iter()
        .map(|e| e.to_vec())
        .collect();
        let mut test = NaivedeBruijnGraph::new(k, b"ACGT");
        for kmer in &kmers {
            assert_eq!(test.insert(kmer), Some(1));
        }
        eprintln!("{:?}", test);
        for kmer in &kmers {
            assert_eq!(test.num(kmer), 1);
        }
        assert_eq!(test.num(b"ACCGTATCGAT"), 0);
        for kmer in &kmers {
            assert_eq!(test.insert(kmer), Some(2));
        }
        for kmer in &kmers {
            assert_eq!(test.num(kmer), 2);
        }
    }
}
