#![feature(test)]
extern crate bio;
extern crate fnv;
extern crate fxhash;
extern crate highway;
extern crate metrohash;
extern crate test;
extern crate twox_hash;

fn main(){
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::hash::Hasher;
    use test::Bencher;
    #[test]
    fn test() {}
    #[bench]
    fn bench(b: &mut test::Bencher) {
        b.iter(|| 1);
    }
    fn open() -> Vec<Vec<u8>> {
        let file = "/grid/ban-m/test.fq";
        let data: Vec<_> = fastq::Reader::from_file(file)
            .unwrap()
            .records()
            .filter_map(|e| e.ok())
            .map(|e| e.seq().to_vec())
            .collect();
        eprintln!("Datasize:{}", data.iter().map(|e| e.len()).sum::<usize>());
        data
    }
    #[bench]
    fn two_hash_16(b: &mut Bencher) {
        let data = open();
        b.iter(|| {
            let res: u32 = data
                .iter()
                .flat_map(|e| e.windows(16))
                .map(|kmer| {
                    let mut hs = twox_hash::XxHash32::with_seed(0);
                    hs.write(kmer);
                    hs.finish() as u32
                })
                .sum();
            test::black_box(res);
        });
    }
    #[bench]
    fn fxhash_16(b: &mut Bencher) {
        let data = open();
        b.iter(|| {
            let res: u32 = data
                .iter()
                .flat_map(|e| e.windows(16))
                .map(fxhash::hash32)
                .sum();
            test::black_box(res);
        })
    }
    #[bench]
    fn fnv_16(b: &mut Bencher) {
        let data = open();
        b.iter(|| {
            let res: u32 = data
                .iter()
                .flat_map(|e| e.windows(16))
                .map(|kmer| {
                    let mut hs = fnv::FnvHasher::default();
                    hs.write(kmer);
                    hs.finish() as u32
                })
                .sum();
            test::black_box(res);
        })
    }
    #[bench]
    fn metrohash_16(b: &mut Bencher) {
        let data = open();
        b.iter(|| {
            let res: u32 = data
                .iter()
                .flat_map(|e| e.windows(16))
                .map(|kmer| {
                    let mut hs = metrohash::MetroHash64::default();
                    hs.write(kmer);
                    hs.finish() as u32
                })
                .sum();
            test::black_box(res);
        })
    }
    #[bench]
    fn highway_16(b: &mut Bencher) {
        let data = open();
        let key = highway::Key([1, 2, 3, 4]);
        use highway::HighwayHash;
        b.iter(|| {
            let res: u32 = data
                .iter()
                .flat_map(|e| e.windows(16))
                .map(|kmer| {
                    let hs = highway::HighwayBuilder::new(&key);
                    hs.hash64(kmer) as u32
                })
                .sum();
            test::black_box(res);
        })
    }
}
