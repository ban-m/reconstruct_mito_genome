pub fn get_start_stop(
    reads: &[crate::EncodedRead],
    contigs: &super::Contigs,
) -> Vec<(u16, Vec<(usize, u32)>)> {
    let mut coverage: Vec<Vec<u32>> = contigs
        .get_last_units()
        .iter()
        .map(|&last| vec![0; last as usize + 1])
        .collect();
    let mut gaps = 0;
    let total = reads.len();
    for read in reads {
        let mut seq = read.seq.iter().skip_while(|e| e.is_gap());
        let head = match seq.next() {
            Some(unit) if unit.is_encode() => unit.encode().unwrap(),
            _ => {
                debug!("Read {} is gappy.", read.id());
                gaps += 1;
                continue;
            }
        };
        coverage[head.contig as usize][head.unit as usize] += 1;
        let mut seq = read.seq.iter().rev().skip_while(|e| e.is_gap());
        let tail = match seq.next() {
            Some(unit) if unit.is_encode() => unit.encode().unwrap(),
            _ => {
                debug!("Read {} is single-unit.", read.id());
                continue;
            }
        };
        coverage[tail.contig as usize][tail.unit as usize] += 1;
    }
    debug!("There are {} gappy reads out of {} reads.", gaps, total);
    coverage
        .into_iter()
        .enumerate()
        .map(|(idx, cov)| {
            let cov: Vec<_> = cov.into_iter().enumerate().collect();
            (idx as u16, cov)
        })
        .collect()
}
