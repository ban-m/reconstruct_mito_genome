extern crate packed_simd;
fn main() {
    use packed_simd::m32x4 as m32s;
    let t1 = m32s::new(true, true, false, false);
    let t2 = m32s::new(true, false, true, false);
    eprintln!("{:?}", t1 & t2);
    eprintln!("{:?}", t1 | t2);
    // use std::time::Instant;
    // let num: Vec<i32> = (0..1600).map(|e| e as i32).collect();
    // let start = Instant::now();
    // for _ in 0..1000 {
    //     let sum = num.iter().fold(0, |x, y| x + y);
    // }
    // eprintln!("Fold:{:?}", Instant::now() - start);
    // let start = Instant::now();
    // for _ in 0..1000 {
    //     let sum = num.iter().sum::<i32>();
    // }
    // eprintln!("Usual:{:?}", Instant::now() - start);
    // let start = Instant::now();
    // use packed_simd::i32x2;
    // for _ in 0..1000 {
    //     let sum = num
    //         .chunks_exact(i32x2::lanes())
    //         .map(|slice| i32x2::from_slice_unaligned(slice))
    //         .fold(i32x2::splat(0), |x, y| x + y)
    //         .wrapping_sum();
    // }
    // eprintln!("i32x2:{:?}", Instant::now() - start);
    // use packed_simd::i32x16;
    // let start = Instant::now();
    // for _ in 0..1000 {
    //     let sum = num
    //         .chunks_exact(i32x16::lanes())
    //         .map(|slice| i32x16::from_slice_unaligned(slice))
    //         .fold(i32x16::splat(0), |x, y| x + y)
    //         .wrapping_sum();
    // }
    // eprintln!("i32x16:{:?}", Instant::now() - start);
}
