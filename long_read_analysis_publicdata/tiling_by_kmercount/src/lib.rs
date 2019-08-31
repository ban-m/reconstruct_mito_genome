extern crate bio;
use bio::data_structures::{bwt, suffix_array::suffix_array};
use std::collections::BTreeMap;
struct KmerDataBase {
    len: usize,
    occ: bwt::Occ,
    bwt: bwt::BWT,
    count_table: [usize; 255],
    kmercount: BTreeMap<usize, usize>,
}

impl KmerDataBase {
    fn new(mut kmer_count: Vec<(Vec<u8>, usize)>) -> Self {
        // Create the 'N'-padding reference sequence.
        Self::radix_sort(&mut kmer_count);
        let (overlapped_seq, kmercount) = Self::overlap_merge(kmer_count);
        let alphabet = bio::alphabets::dna::n_alphabet();
        let sa = suffix_array(&overlapped_seq);
        let bwt = bwt::bwt(&overlapped_seq, &sa);
        let occ = bwt::Occ::new(&bwt, 1, &alphabet);
        let count_table = Self::count_table(&overlapped_seq);
        let len = overlapped_seq.len();
        Self {
            len,
            occ,
            bwt,
            count_table,
            kmercount,
        }
    }
    fn radix_sort(kmer_count:&mut Vec<(Vec<u8>, usize)>) {
        kmer_count.sort_by_key(|e|e.0);
    }
    fn overlap_merge(kmer_count: Vec<(Vec<u8>, usize)>) -> (Vec<u8>, BTreeMap<usize, usize>) {
        // should be end with '$'.
        let complete_sequence = vec![];
        let position_to_count = BTreeMap::new();
        let mut previous_kmer = vec![];
        for (kmer,count) in kmer_count{
            
        }
        (complete_sequence, position_to_count)
    }
    fn count_table(seq: &[u8]) -> [usize; 255] {
        [0; 255]
    }
}
