extern crate bio;
extern crate rayon;
use bio::io::fasta;
use rayon::prelude::*;

fn to_index(kmer: &[u8], k: usize) -> usize {
    kmer.iter()
        .enumerate()
        .filter_map(|(dig, e)| match e {
            b'a' | b'A' => Some((k - dig, 0)),
            b'c' | b'C' => Some((k - dig, 1)),
            b'g' | b'G' => Some((k - dig, 2)),
            b't' | b'T' => Some((k - dig, 3)),
            _ => None,
        })
        .map(|(dig, num)| num * (4_usize.pow(dig as u32 - 1)))
        .sum()
}

#[test]
fn test_to_index() {
    let kmer = b"AAAAA";
    assert_eq!(to_index(kmer, 5), 0);
    let kmer = b"AAAAC";
    assert_eq!(to_index(kmer, 5), 1);
    let kmer = b"AAAAT";
    assert_eq!(to_index(kmer, 5), 3);
    let kmer = b"TAAAA";
    assert_eq!(to_index(kmer, 5), 4usize.pow(4) * 3);
    let kmer = b"TAAAC";
    assert_eq!(to_index(kmer, 5), 4usize.pow(4) * 3 + 1);
}

#[test]
fn test_into_kmer() {
    let k: usize = 5;
    assert_eq!(into_kmer(0, k), "AAAAA");
    assert_eq!(into_kmer(1, k), "AAAAC");
    assert_eq!(into_kmer(3, k), "AAAAT");
    assert_eq!(into_kmer(4usize.pow(4) * 3, k), "TAAAA");
    let k: usize = 10;
    for e in 0..4usize.pow(k as u32) {
        let index = to_index(into_kmer(e, k).as_bytes(), k);
        assert_eq!(e, index);
    }
}

fn convert_to_vector_freq(record: &fasta::Record, k: usize) -> Vec<f64> {
    let len = record.seq().len() as f64;
    let mut freq = vec![0u16; 4_usize.pow(k as u32)];
    for kmer in record.seq().windows(k) {
        freq[to_index(kmer, k)] += 1;
    }
    freq.into_iter().map(|e| e as f64 / len).collect()
}

fn into_kmer(e: usize, k: usize) -> String {
    eprintln!("{:b}\t{}", e, k);
    let mut kmer = String::new();
    for i in 0..k {
        let c = e >> 2 * (k - i - 1) & 0b11;
        eprintln!("{:b} After {} bit shift", c, k - i - 1);
        match c {
            0 => kmer.push('A'),
            1 => kmer.push('C'),
            2 => kmer.push('G'),
            3 => kmer.push('T'),
            _ => {}
        };
    }
    kmer
}

fn all_kmers(k: usize) -> Vec<String> {
    (0..(4_usize.pow(k as u32)))
        .map(|e| into_kmer(e, k))
        .collect()
}

fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let k: usize = args[2].parse().unwrap();
    let records: Vec<_> = fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
        .records()
        .filter_map(|e| e.ok())
        .collect();

    let freq_vectors: Vec<_> = records
        .into_par_iter()
        .map(|rec| (rec.id().to_string(), convert_to_vector_freq(&rec, k)))
        .collect();
    print!("id\t{}", all_kmers(k).join("\t"));
    for (id, freq) in freq_vectors {
        let freq:Vec<_> = freq.into_iter().map(|e|format!("{}",e)).collect();
        let freq = freq.join("\t");
        println!("{}\t{}",id,freq);
    }
    Ok(())
}
