use last_tiling::unit::*;
use last_tiling::Contigs;
use std::collections::HashSet;

const READ_NUM: usize = 8;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum CriticalRegion {
    BR(BroadRepeat),
    CR(ConfluentRegion),
    JP(JointPoint),
}

/// A trait to classify a given read.
pub trait ReadClassify {
    /// Whether the read spans `self`.
    fn is_spanned_by(&self, r: &EncodedRead) -> bool;
}

#[derive(Debug, Clone)]
pub struct BroadRepeat {
    contig1: Position,
    contig2: Position,
}

/// Return the maximum/minumum index of encoded read's unit with contig of `contig`.
pub fn get_max_min_unit(r: &EncodedRead, contig: u16) -> (u16, u16) {
    r.seq()
        .iter()
        .filter_map(|e| match e {
            ChunkedUnit::En(e) if e.contig == contig => Some(e.unit),
            _ => None,
        })
        .fold((std::u16::MAX, std::u16::MIN), |(max, min), x| {
            if x < min {
                (x, max)
            } else if max < x {
                (min, x)
            } else {
                (min, max)
            }
        })
}

impl ReadClassify for BroadRepeat {
    // The broad repeat is like the below.
    // -----A||||||B----
    // -----C||||||D----
    // There can be *6* pattern of spanning
    // Let them 1:A<->B, 2:A<->D, 3:A<->C, 4:C<->D, 5:C<->B, 6:D<->B).
    fn is_spanned_by(&self, r: &EncodedRead) -> bool {
        let c1 = self.contig1.contig;
        let c2 = self.contig2.contig;
        let (c1_min, c1_max) = get_max_min_unit(r, c1);
        let (c2_min, c2_max) = get_max_min_unit(r, c2);
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        let span1 = c1_min < c1_s && c1_e < c1_max;
        let span2 = c1_min < c1_s && c2_e < c2_max;
        let span3 = c1_min < c1_s && c2_min < c2_s;
        let span4 = c2_min < c2_s && c2_e < c2_max;
        let span5 = c2_min < c2_s && c1_e < c1_max;
        let span6 = c2_e < c2_max && c1_e < c1_max;
        span1 || span2 || span3 || span4 || span5 || span6
    }
}

impl BroadRepeat {
    fn new(contig1: Position, contig2: Position) -> Self {
        Self { contig1, contig2 }
    }
}

#[derive(Debug, Clone)]
pub struct ConfluentRegion {
    start_1: Position,
    start_2: Position,
    end_1: Position,
    end_2: Position,
}

#[derive(Debug, Clone)]
pub struct JointPoint {
    start_1: Position,
    start_2: Position,
    end_1: Position,
    end_2: Position,
}

/// The position at contigs.
#[derive(Debug, Clone)]
pub struct Position {
    contig: u16,
    start_unit: u16,
    end_unit: u16,
}

impl Position {
    fn new(contig: u16, s: u16, e: u16) -> Self {
        let (start_unit, end_unit) = (s.min(e), s.max(e));
        Self {
            contig,
            start_unit,
            end_unit,
        }
    }
    fn range(&self) -> (u16, u16) {
        (self.start_unit, self.end_unit)
    }
}

/// Return critical regions.
pub fn critical_regions(reads: &[EncodedRead], contigs: &Contigs) -> Vec<BroadRepeat> {
    let num_of_contig = contigs.get_num_of_contigs() as u16;
    let mut regions = vec![];
    for from in 0..num_of_contig {
        if let Some(res) = critical_region_within(from, reads, contigs) {
            regions.extend(res);
        }
        for to in from + 1..num_of_contig {
            if let Some(res) = critical_regions_between(from, to, reads, contigs) {
                regions.extend(res);
            }
        }
    }
    regions
}

fn critical_region_within(
    contig: u16,
    reads: &[EncodedRead],
    contigs: &Contigs,
) -> Option<Vec<BroadRepeat>> {
    let last_unit = contigs.get_last_unit(contig)? as usize;
    let mut inner_count: Vec<Vec<&EncodedRead>> = vec![vec![]; last_unit + 1];
    for read in reads {
        // head seeking
        let mut units = read.seq().iter().skip_while(|e| e.is_gap());
        let mut prev = match units.next()? {
            ChunkedUnit::En(e) => e,
            _ => unreachable!(),
        };
        let mut num_gap = 0;
        let hg = last_tiling::UNIT_SIZE / 2;
        for unit in units {
            match unit {
                ChunkedUnit::En(e) if e.contig == contig => {
                    let is_forward =
                        e.unit == prev.unit + num_gap + 1 || e.unit + num_gap == prev.unit + 1;
                    let is_reverse =
                        e.unit + 1 == prev.unit + num_gap || e.unit + 1 + num_gap == prev.unit;
                    num_gap = 0;
                    if !is_forward && !is_reverse {
                        // Gapping. Put flag.
                        inner_count[e.unit as usize].push(read);
                        inner_count[prev.unit as usize].push(read);
                    }
                    prev = e;
                }
                ChunkedUnit::Gap(g) => num_gap += ((g.len() + hg) / last_tiling::UNIT_SIZE) as u16,
                _ => num_gap += 1,
            }
        }
    }
    let inner_region = calculate_average_more_than(inner_count.iter().map(|e| e.len()).collect());
    let mut result = vec![];
    for i in 0..inner_region.len() {
        for j in i..inner_region.len() {
            let &(s, t) = &inner_region[i];
            let &(x, y) = &inner_region[j];
            let c1 = Position::new(contig, s, t);
            let c2 = Position::new(contig, x, y);
            let f: HashSet<_> = inner_count[s as usize..t as usize].iter().collect();
            let g: HashSet<_> = inner_count[x as usize..y as usize].iter().collect();
            if f.intersection(&g).count() > READ_NUM {
                result.push(BroadRepeat::new(c1, c2))
            }
        }
    }
    Some(result)
}

fn critical_regions_between(
    from: u16,
    to: u16,
    reads: &[EncodedRead],
    contigs: &Contigs,
) -> Option<Vec<BroadRepeat>> {
    // Enumerate critical regions from `from` contig to `to` contig.
    let from_last_unit = contigs.get_last_unit(from)?;
    let to_last_unit = contigs.get_last_unit(to)?;
    let mut from_count: Vec<Vec<&EncodedRead>> = vec![vec![]; from_last_unit as usize + 1];
    let mut to_count: Vec<Vec<&EncodedRead>> = vec![vec![]; to_last_unit as usize + 1];
    for read in reads {
        let mut units = read.seq().iter().filter_map(|e| match e {
            ChunkedUnit::En(ref encode) if encode.contig == from || encode.contig == to => {
                Some(encode)
            }
            _ => None,
        });
        let mut prev = units.next()?;
        for unit in units {
            if prev.contig == from && unit.contig == to {
                from_count[prev.unit as usize].push(read);
                to_count[unit.unit as usize].push(read);
            } else if prev.contig == to && unit.contig == from {
                from_count[unit.unit as usize].push(read);
                to_count[prev.unit as usize].push(read);
            }
            prev = unit;
        }
    }
    let from_region = calculate_average_more_than(from_count.iter().map(|e| e.len()).collect());
    let to_region = calculate_average_more_than(to_count.iter().map(|e| e.len()).collect());
    let mut regions = vec![];
    for &(s, t) in &from_region {
        for &(x, y) in &to_region {
            let c1 = Position::new(from, s, t);
            let c2 = Position::new(to, x, y);
            // Determine whether these pairs would be a paired region.
            let f: HashSet<&EncodedRead> = from_count[s as usize..t as usize]
                .iter()
                .flat_map(|e| e.iter())
                .map(|&e| e)
                .collect();
            let g: HashSet<&EncodedRead> = to_count[x as usize..y as usize]
                .iter()
                .flat_map(|e| e.iter())
                .map(|&e| e)
                .collect();
            if g.intersection(&f).count() > READ_NUM {
                regions.push(BroadRepeat::new(c1, c2));
                // if 0 < s && t < from_last_unit + 1 && 0 < x && y < to_last_unit + 1 {
                //     regions.push()
                // } else if (s == 0 || t == from_last_unit + 1) && (0 < x && y < to_last_unit + 1) {
                //     // Confluent Region
                // } else if (0 < s && t < from_last_unit + 1) && (x == 0 || y == to_last_unit + 1) {
                //     // Confluent Region
                // } else {
                //     // Joit Point
                // }
            }
        }
    }
    Some(regions)
}

fn calculate_average_more_than(input: Vec<usize>) -> Vec<(u16, u16)> {
    // Calculate maximal region with average more than READ_NUM.
    // 1. Extract average.
    let input: Vec<i32> = input
        .into_iter()
        .map(|e| e as i32 - READ_NUM as i32)
        .collect();
    // 2. Cumsum. If c[i] < c[j], [i,j) would satisfy the condition.
    // Here, c[j] is the sum up to index j-1 !! Thus,  always c[0] == 0.
    // And, c[input.len()] is valid. This is because input[0..input.len()] should be valid.
    let (_total, cumsum): (i32, Vec<i32>) =
        input.iter().fold((0, vec![0]), |(acc, mut cumsum), &x| {
            cumsum.push(acc + x);
            (x + acc, cumsum)
        });
    assert!(cumsum.len() == input.len() + 1);
    // rmin[i] = min_{j=0}^i c[j].
    // lmax[i] = max_{j=i}^n c[j].
    let mut rmin = vec![cumsum[0]; cumsum.len()];
    for i in 1..cumsum.len() {
        rmin[i] = rmin[i - 1].min(cumsum[i]);
    }
    let mut lmax = vec![*cumsum.last().unwrap(); cumsum.len()];
    for i in (0..cumsum.len() - 1).rev() {
        lmax[i] = lmax[i + 1].max(cumsum[i]);
    }
    let (mut i, mut j) = (0, 0);
    let mut result = vec![];
    while i <= input.len() && j <= input.len() {
        if rmin[i] < lmax[j] {
            j += 1;
        } else {
            // Output j-1 and i.
            if 0 < j {
                result.push((i as u16, j as u16));
            }
            while lmax[j] <= rmin[i] {
                i += 1;
            }
        }
    }
    result
}
