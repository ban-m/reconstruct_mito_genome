extern crate bio_utils;
use bio_utils::maf;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let records: Vec<_> = maf::Reader::from_file(&args[1])?
        .records()
        .filter_map(|e| e.ok())
        .filter(|e| match e.score() {
            Some(score) => score > 500.,
            None => false,
        })
        .map(|e| {
            let seq = e.sequence();
            assert!(seq.len() >= 2);
            let (sum, mismat, del) = seq[0].text().iter().zip(seq[1].text().iter()).fold(
                (0, 0, 0),
                |(s, m, d), (a, b)| {
                    if a == &b'-' || b == &b'-' {
                        (s + 1, m, d + 1)
                    } else if a.to_ascii_uppercase() != b.to_ascii_uppercase() {
                        (s + 1, m + 1, d)
                    } else {
                        (s + 1, m, d)
                    }
                },
            );
            (mismat as f64 / sum as f64, del as f64 / sum as f64)
        })
        .collect();
    let len = records.len() as f64;
    let (mismat_ave, del_ave) = records
        .iter()
        .fold((0., 0.), |(x, y), (m, d)| (x + m, y + d));
    let mismat_ave = mismat_ave / len;
    let del_ave = del_ave / len;
    let records: Vec<_> = records.into_iter().map(|(d, m)| d + m).collect();
    let max = records
        .iter()
        .fold(std::f64::MIN, |max, &x| if max < x { x } else { max });
    let min = records
        .iter()
        .fold(std::f64::MAX, |min, &x| if min < x { min } else { x });
    let ave = records.iter().fold(0., |s, &x| s + x) / records.len() as f64;
    println!("{:?}", records);
    println!("Max\tMin\tAve\tdel\tmismat");
    println!("{},{},{},{},{}",max, min, ave, del_ave, mismat_ave);
    Ok(())
}
