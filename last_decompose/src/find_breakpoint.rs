use super::UNIT_SIZE;
use last_tiling::unit::*;
use last_tiling::Contigs;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
/// The threshould used to determine whether two regions would be regarded as pairs.
const READ_NUM: usize = 4;
/// How many reads are needed to construct a separate class.
/// Note that it should be more than 1, since there is a
/// danger of chimeric reads.
const CLASS_THR: usize = 4;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CriticalRegion {
    CP(ContigPair),
    // CR(ConfluentRegion),
    RJ(RepeatJunction),
}

// A trait to classify a given read.
// (read class, read index, read itself)
type ReadWithClassAndIndex<'a> = (usize, usize, &'a EncodedRead);
pub trait ReadClassify {
    /// Whether the read spans `self`.
    fn is_spanned_by(&self, r: &EncodedRead) -> bool;
    /// Whether the read touches `self`. Return false if the read is contained.
    fn overlaps_with(&self, r: &EncodedRead) -> bool;
    /// Wether `self` contains the read.
    fn contains(&self, r: &EncodedRead) -> bool;
    /// Separate (spanned) reads into corresponding slots.
    /// It automatically remove some of the reads supporting very weak connections.
    /// This method should not invoke heavy procs such as HMMs or other prediction methods.
    /// For example, `self` should use the most naive approach such as just distributing reads.
    fn separate_reads_into_clusters<'a>(
        &self,
        reads: Vec<(usize, &'a EncodedRead)>,
    ) -> (usize, Vec<ReadWithClassAndIndex<'a>>);
    /// The shortest distance from r to self.
    fn distance(&self, r: &EncodedRead) -> usize;
}

impl ReadClassify for CriticalRegion {
    fn overlaps_with(&self, r: &EncodedRead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.overlaps_with(r),
            CriticalRegion::RJ(ref rj) => rj.overlaps_with(r),
        }
    }
    fn is_spanned_by(&self, r: &EncodedRead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.is_spanned_by(r),
            CriticalRegion::RJ(ref rj) => rj.is_spanned_by(r),
        }
    }
    fn contains(&self, r: &EncodedRead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.contains(r),
            CriticalRegion::RJ(ref rj) => rj.contains(r),
        }
    }
    fn separate_reads_into_clusters<'a>(
        &self,
        reads: Vec<(usize, &'a EncodedRead)>,
    ) -> (usize, Vec<ReadWithClassAndIndex<'a>>) {
        match self {
            CriticalRegion::CP(ref cp) => cp.separate_reads_into_clusters(reads),
            CriticalRegion::RJ(ref rj) => rj.separate_reads_into_clusters(reads),
        }
    }
    fn distance(&self, r: &EncodedRead) -> usize {
        match self {
            CriticalRegion::CP(ref cp) => cp.distance(r),
            CriticalRegion::RJ(ref rj) => rj.distance(r),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigPair {
    contig1: Position,
    contig2: Position,
}

/// Return the maximum/minumum index of encoded read's unit with contig of `contig`.
/// Note that if the start unit is 0-th unit of `contig`, followed by large gap>UNIT_SIZE
/// It assigns -1.
/// And if the reads ended at i-th unit of `contig`, followed by large gap>UNIT_SIZE,
/// it assigns i+1.
pub fn get_max_min_unit(r: &EncodedRead, contig: u16) -> (i32, i32) {
    let (mut min, mut max) = (std::i32::MAX, std::i32::MIN);
    let mut last_met_unitnum = None;
    let mut had_not_met_encode = true;
    let mut last_unit_was_large_gap = false;
    let mut first_unit_was_large_gap = false;
    for unit in r.seq() {
        match unit {
            ChunkedUnit::En(e) if e.contig == contig => {
                if had_not_met_encode && first_unit_was_large_gap {
                    min = e.unit as i32 - 1;
                    max = e.unit as i32 + 1;
                } else {
                    min = min.min(e.unit as i32);
                    max = max.max(e.unit as i32);
                }
                had_not_met_encode = false;
                last_unit_was_large_gap = false;
                last_met_unitnum = Some(e.unit as i32);
            }
            ChunkedUnit::En(_) => {
                had_not_met_encode = false;
                last_unit_was_large_gap = false;
            }
            ChunkedUnit::Gap(g) if g.len() > UNIT_SIZE => {
                last_unit_was_large_gap = true;
                first_unit_was_large_gap = true;
            }
            _ => {}
        }
    }
    // It indicate that last encoded unit is the wanted one, and followed by long gap.
    // Increment the maximum number of unit by 1.
    if last_met_unitnum == Some(max) && last_unit_was_large_gap {
        max += 1;
    } else if last_met_unitnum == Some(min) && last_unit_was_large_gap {
        min -= 1;
    }
    (min, max)
}

// A helper method for ReadClassify for ContigPair.
// This function takes the classed reads, assesses if
// there is enough number of reads.
// If it has, this function push these reads into result bucket.
fn extend_and_increment<'a>(
    result: &mut Vec<ReadWithClassAndIndex<'a>>,
    class: usize,
    reads: Vec<&(usize, &'a EncodedRead)>,
) -> usize {
    if reads.len() > CLASS_THR {
        result.extend(reads.into_iter().map(|&(idx, r)| (class, idx, r)));
        class + 1
    } else {
        class
    }
}

impl ReadClassify for ContigPair {
    // The broad repeat is like the below.
    // -----A||||||B----
    // -----C||||||D----
    // There can be *6* pattern of spanning
    // Let them 1:A<->B, 2:A<->D, 3:A<->C, 4:C<->D, 5:C<->B, 6:D<->B).
    // Even thought Either A<->C or A<->D is feasible, since the
    // strand information, I check the all possibilities
    // such as non-stringient recombinations.
    // Remark: It can handle the tail gap correctly.
    // Formally speaking, it converts all the variables into signed variable,
    // and assign -1 to head gap, +1 to tail gap.
    fn is_spanned_by(&self, r: &EncodedRead) -> bool {
        let c1 = self.contig1.contig;
        let c2 = self.contig2.contig;
        let (c1_min, c1_max) = get_max_min_unit(r, c1);
        let (c2_min, c2_max) = get_max_min_unit(r, c2);
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        let span1 = (c1_min < c1_s || c2_min < c2_s) && (c1_e < c1_max || c2_e < c2_max);
        let span2 = c1_min < c1_s && c2_min < c2_s;
        let span3 = c2_e < c2_max && c1_e < c1_max;
        span1 || span2 || span3
    }
    // Overlapping can be detected similary.
    // Note that overlapping here means that
    // the read "crosses" boundary or not.
    // Shortly, if the read spans this region, it overlaps.
    // Thus, it is sufficient and nessesary condition to
    // check whther min < boundary < max holds for one of four
    // breakpoints.
    fn overlaps_with(&self, r: &EncodedRead) -> bool {
        let (c1_min, c1_max) = get_max_min_unit(r, self.contig1.contig);
        let (c2_min, c2_max) = get_max_min_unit(r, self.contig2.contig);
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        // The code below is correst, as c*_min would be std::i32::MAX
        // if there is no contig * unit in the read `r`, and c*_max would be std::i32::MIN
        // if there is no contig * unit.
        c1_min <= c1_s && c1_s <= c1_max
            || c1_min <= c1_e && c1_e <= c1_max
            || c2_min <= c2_s && c2_s <= c2_max
            || c2_min <= c2_e && c2_e <= c2_max
    }
    // Assign each reads into its color. Note that the
    // color should be at most four, depending of the spanning patterns.
    // Returns the number of cluster, with a vector of (class,index,read) tuples.
    fn separate_reads_into_clusters<'a>(
        &self,
        reads: Vec<(usize, &'a EncodedRead)>,
    ) -> (usize, Vec<ReadWithClassAndIndex<'a>>) {
        // Determine class.
        assert!(reads.iter().all(|&(_idx, r)| self.is_spanned_by(r)));
        let c1 = self.contig1.contig;
        let c2 = self.contig2.contig;
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        // TODO: refactoring this function. Currently very agry.
        let mut result = vec![];
        let mut cluster_num = 0;
        {
            let c1s_to_c1e: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (c1_min, c1_max) = get_max_min_unit(r, c1);
                    c1_min < c1_s && c1_e < c1_max
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c1s_to_c1e);
        }
        {
            let c2s_to_c2e: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (c2_min, c2_max) = get_max_min_unit(r, c2);
                    c2_min < c2_s && c2_e < c2_max
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c2s_to_c2e);
        }
        {
            let c1s_to_c2e: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (c1_min, _) = get_max_min_unit(r, c1);
                    let (_, c2_max) = get_max_min_unit(r, c2);
                    c1_min < c1_s && c2_e < c2_max
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c1s_to_c2e);
        }
        {
            let c2s_to_c1e: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (_, c1_max) = get_max_min_unit(r, c1);
                    let (c2_min, _) = get_max_min_unit(r, c2);
                    c2_min < c2_e && c1_e < c1_max
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c2s_to_c1e);
        }
        {
            let c1s_to_c2s: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (c1_min, _) = get_max_min_unit(r, c1);
                    let (c2_min, _) = get_max_min_unit(r, c2);
                    c1_min < c1_s && c2_min < c2_s
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c1s_to_c2s);
        }
        {
            let c1e_to_c2e: Vec<_> = reads
                .iter()
                .filter(|&(_, r)| {
                    let (_, c1_max) = get_max_min_unit(r, c1);
                    let (_, c2_max) = get_max_min_unit(r, c2);
                    c1_e < c1_max && c2_e < c2_max
                })
                .collect();
            cluster_num = extend_and_increment(&mut result, cluster_num, c1e_to_c2e);
        }
        (cluster_num, result)
    }
    // Return whether self contains the read.
    // It is enough to check s < start && e < end holds
    // for contig 1 or contig 2.
    // Remember that, if there's no unit,
    // get_max_min_unit(,) would return (MAX,MIN), where
    // the condition above never holds.
    fn contains(&self, r: &EncodedRead) -> bool {
        let (c1_s, c1_e) = self.contig1.range();
        let (c1_min, c1_max) = get_max_min_unit(r, self.contig1.contig);
        let contained_in_c1 = c1_s < c1_min && c1_max < c1_e;
        let (c2_s, c2_e) = self.contig1.range();
        let (c2_min, c2_max) = get_max_min_unit(r, self.contig2.contig);
        let contained_in_c2 = c2_s < c2_min && c2_max < c2_e;
        contained_in_c1 || contained_in_c2
    }
    // Calculate the distance between the read and self.
    // It can be assumed that the read is not contained nor not overlapping with self.
    // For the reasen, search for the location where this function is called.
    fn distance(&self, r: &EncodedRead) -> usize {
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        let (c1_min, c1_max) = get_max_min_unit(r, self.contig1.contig);
        let (c2_min, c2_max) = get_max_min_unit(r, self.contig2.contig);
        [c1_s - c1_max, c1_e - c1_min, c2_s - c2_max, c2_e - c2_min]
            .into_iter()
            .copied()
            .map(i32::abs)
            .min()
            .unwrap() as usize
    }
}

impl ContigPair {
    pub fn new(contig1: Position, contig2: Position) -> Self {
        Self { contig1, contig2 }
    }
    pub fn get_competitive_pair(&self, ovlp: &[(usize, &EncodedRead)]) -> Vec<[usize; 2]> {
        // C1:-----1=======2------
        // C2:-----3=======4------
        // (A,C) vs (B,D) or (A,D) vs (B,C)
        // It can be determined by counting concording reads(forward in both C1 and C2, or reverse in both C1 and C2)
        // and non-concording reads(forward in C1, reverse in C2 and vice versa).
        let (concord, disconcord) = ovlp
            .iter()
            .filter_map(|(_, r)| {
                Some((
                    r.is_forward_wrt(self.contig1.contig)?,
                    r.is_forward_wrt(self.contig2.contig)?,
                ))
            })
            .fold(
                (0, 0),
                |(c, d), (c1, c2)| if c1 == c2 { (c + 1, d) } else { (c, d + 1) },
            );
        if concord > disconcord {
            vec![[1, 3], [2, 4]]
        } else {
            vec![[1, 4], [2, 3]]
        }
    }
    // TODO
    pub fn clustering_reads_into_points<'a>(
        &self,
        reads: &[(usize, &'a EncodedRead)],
    ) -> (usize, Vec<Vec<(usize, &'a EncodedRead)>>) {
        // Distribute reads into four points listed below:
        // C1:-----1=======2------
        // C2:-----3=======4------
        let num_of_cluster = 4;
        let mut result = vec![vec![]; 4];
        let (c1_s, c1_e) = self.contig1.range();
        let (c2_s, c2_e) = self.contig2.range();
        for &(idx, read) in reads {
            // Determine point.
            let (c1min, c1max) = get_max_min_unit(read, self.contig1.contig);
            let (c2min, c2max) = get_max_min_unit(read, self.contig2.contig);
            if c1min < c1_s && c1_s < c1max {
                result[0].push((idx, read));
            } else if c1min < c1_e && c1_e < c1max {
                result[1].push((idx, read));
            } else if c2min < c2_s && c2_s < c2max {
                result[2].push((idx, read));
            } else if c2min < c2_e && c2_e < c2max {
                result[3].push((idx, read));
            }
        }
        (num_of_cluster, result)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepeatJunction {
    // Contains start and end position.
    pos: Position,
}

impl RepeatJunction {
    fn new(contig: u16, start: u16, end: u16) -> Self {
        let pos = Position::new(contig, start, end);
        Self { pos }
    }
}

impl ReadClassify for RepeatJunction {
    fn is_spanned_by(&self, r: &EncodedRead) -> bool {
        let (min, max) = get_max_min_unit(r, self.pos.contig);
        let (start, end) = self.pos.range();
        // Sometimes the start position and end position is
        // exactly the boundary, thus, it is indeed needed.
        let some_unit = r
            .seq()
            .iter()
            .filter_map(|e| e.encode())
            .filter(|e| e.contig != self.pos.contig)
            .count()
            != 0;
        min <= start && end <= max && some_unit
    }
    fn overlaps_with(&self, r: &EncodedRead) -> bool {
        let some_unit = r
            .seq()
            .iter()
            .filter_map(|e| e.encode())
            .filter(|e| e.contig != self.pos.contig)
            .count()
            != 0;
        let (min, max) = get_max_min_unit(r, self.pos.contig);
        min != std::i32::MAX && max != std::i32::MIN && some_unit
    }
    fn contains(&self, r: &EncodedRead) -> bool {
        let no_other_unit = r
            .seq()
            .iter()
            .filter_map(|e| e.encode())
            .filter(|e| e.contig != self.pos.contig)
            .count()
            == 0;
        let (min,max) = get_max_min_unit(r,self.pos.contig);
        let (start,end) = self.pos.range();
        start  < min && max < end && no_other_unit
    }
    fn separate_reads_into_clusters<'a>(
        &self,
        _reads: Vec<(usize, &'a EncodedRead)>,
    ) -> (usize, Vec<ReadWithClassAndIndex<'a>>) {
        (0, vec![])
    }
    // TODO
    fn distance(&self, _r: &EncodedRead) -> usize {
        0
    }
}

/// The position at contigs.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    fn range(&self) -> (i32, i32) {
        (self.start_unit as i32, self.end_unit as i32)
    }
}

/// Return critical regions.
pub fn critical_regions(reads: &[EncodedRead], contigs: &Contigs) -> Vec<CriticalRegion> {
    let num_of_contig = contigs.get_num_of_contigs() as u16;
    let mut regions = vec![];
    debug!("There are {} contigs", num_of_contig);
    for c in contigs.names().iter() {
        debug!("->{}({}len)", c, contigs.get(c).unwrap().len());
    }
    for from in 0..num_of_contig {
        if let Some(res) = critical_region_within(from, reads, contigs) {
            for c in &res {
                debug!("{:?}", c);
            }
            regions.extend(res);
        }
        debug!("Inner end");
        for to in from + 1..num_of_contig {
            if let Some(res) = critical_regions_between(from, to, reads, contigs) {
                for c in &res {
                    debug!("{:?}", c);
                }
                regions.extend(res);
            }
        }
        debug!("Pair end");
        if let Some(res) = critical_regions_repeat_junction(from, reads, contigs) {
            for c in &res {
                debug!("{:?}", c);
            }
            regions.extend(res);
        }
    }
    regions
}

fn critical_region_within(
    contig: u16,
    reads: &[EncodedRead],
    contigs: &Contigs,
) -> Option<Vec<CriticalRegion>> {
    if contigs.is_repeat(contig) {
        return None;
    }
    debug!("{}-th contig, self critical region.", contig);
    let last_unit = contigs.get_last_unit(contig)? as usize;
    let mut inner_count: Vec<Vec<&EncodedRead>> = vec![vec![]; last_unit + 1];
    for read in reads.iter().filter(|r| r.has(contig)) {
        // head seeking
        let mut units = read.seq().iter().skip_while(|e| e.is_gap());
        let mut prev = match units.next()? {
            ChunkedUnit::En(e) => e,
            _ => unreachable!(),
        };
        let mut num_gap = 0;
        for unit in units {
            match unit {
                ChunkedUnit::En(e) if e.contig == contig => {
                    let is_forward =
                        e.unit == prev.unit + num_gap + 1 || e.unit + num_gap == prev.unit + 1;
                    let is_reverse =
                        e.unit + 1 == prev.unit + num_gap || e.unit + 1 + num_gap == prev.unit;
                    // if prev.contig == contig && (prev.unit + 1 != e.unit && prev.unit != e.unit + 1)
                    // {
                    //     debug!("{}-{}(gap:{}),{}", prev, e, num_gap, !is_forward && !is_reverse);
                    // }
                    if !is_forward && !is_reverse && prev.contig == contig {
                        // debug!("{}-{}(gap:{})", prev, e, num_gap);
                        inner_count[e.unit as usize].push(read);
                        inner_count[prev.unit as usize].push(read);
                    }
                    num_gap = 0;
                    prev = e;
                }
                ChunkedUnit::Gap(g) => num_gap += (g.len() / last_tiling::UNIT_SIZE) as u16,
                _ => num_gap += 1,
            }
        }
    }
    debug!("Profiled!({}len)\nPosition\tCount", inner_count.len());
    // for (idx, count) in inner_count
    //     .iter()
    //     .map(|e| e.len())
    //     .enumerate()
    //     .filter(|&(_, len)| len != 0)
    // {
    //     debug!("{}\t{}", idx, count);
    // }
    let inner_counts = inner_count.iter().map(|e| e.len()).collect::<Vec<_>>();
    let inner_region = calculate_average_more_than(&inner_counts, READ_NUM as i32);
    let inner_region = merge_overlap_and_remove_contained(inner_region);
    let inner_region: Vec<_> = inner_region
        .into_iter()
        .map(|(s, e)| tighten_up(s, e, &inner_counts))
        .map(|(s, e)| (s as u16, e as u16))
        .collect();
    debug!("Aggregated!");
    debug!("Start\tEnd");
    // for &(s, t) in inner_region.iter() {
    //     debug!("{}\t{}", s, t);
    // }
    let mut result = vec![];
    for i in 0..inner_region.len() {
        for j in (i + 1)..inner_region.len() {
            let &(s, t) = &inner_region[i];
            let &(x, y) = &inner_region[j];
            if (s <= x && y <= t) || (x <= s && t <= y) {
                continue;
            }
            let c1 = Position::new(contig, s, t);
            let c2 = Position::new(contig, x, y);
            let f: HashSet<_> = inner_count[s as usize..t as usize]
                .iter()
                .flat_map(|e| e.iter())
                .collect();
            let g: HashSet<_> = inner_count[x as usize..y as usize]
                .iter()
                .flat_map(|e| e.iter())
                .collect();
            let in_common = f.intersection(&g).count();
            // debug!(
            //     "{:?}({})-{:?}({}), {} in common.",
            //     (s, t),
            //     f.len(),
            //     (x, y),
            //     g.len(),
            //     in_common
            // );
            if in_common > READ_NUM {
                result.push(CriticalRegion::CP(ContigPair::new(c1, c2)))
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
) -> Option<Vec<CriticalRegion>> {
    // Enumerate critical regions from `from` contig to `to` contig.
    if contigs.is_repeat(from) || contigs.is_repeat(to) {
        debug!("Either {} or {} are repetitive", from, to);
        return None;
    }
    let from_last_unit = contigs.get_last_unit(from)?;
    let to_last_unit = contigs.get_last_unit(to)?;
    let mut from_count: Vec<Vec<&EncodedRead>> = vec![vec![]; from_last_unit as usize + 1];
    let mut to_count: Vec<Vec<&EncodedRead>> = vec![vec![]; to_last_unit as usize + 1];
    for read in reads.iter().filter(|r| r.has(from) && r.has(to)) {
        let mut units = read.seq().iter().filter_map(|e| match e {
            ChunkedUnit::En(ref encode) if encode.contig == from || encode.contig == to => {
                Some(encode)
            }
            _ => None,
        });
        let mut prev = units.next().unwrap();
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
    debug!("btw {}-{}", from, to);
    let from_region: Vec<_> = {
        let raw_count: Vec<_> = from_count.iter().map(|e| e.len()).collect();
        debug!("Row count of {}", from);
        for (idx, count) in raw_count
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 0)
        {
            debug!("{}\t{}", idx, count);
        }
        let res = calculate_average_more_than(&raw_count, READ_NUM as i32);
        for &(s, e) in &res {
            debug!("FROM:[{},{})", s, e);
        }
        let res = merge_overlap_and_remove_contained(res);
        res.into_iter()
            .map(|(s, e)| tighten_up(s, e, &raw_count))
            .map(|(s, e)| (s as u16, e as u16))
            .collect()
    };
    let to_region: Vec<_> = {
        let raw_count: Vec<_> = to_count.iter().map(|e| e.len()).collect();
        debug!("Row count of {}", to);
        for (idx, count) in raw_count
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 0)
        {
            debug!("{}\t{}", idx, count);
        }
        let res = calculate_average_more_than(&raw_count, READ_NUM as i32);
        for &(s, e) in &res {
            debug!("TO:[{},{})", s, e);
        }
        let res = merge_overlap_and_remove_contained(res);
        res.into_iter()
            .map(|(s, e)| tighten_up(s, e, &raw_count))
            .map(|(s, e)| (s as u16, e as u16))
            .collect()
    };
    debug!("{:?}", from_region);
    debug!("{:?}", to_region);
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
                regions.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
        }
    }
    Some(regions)
}

fn critical_regions_repeat_junction(
    from: u16,
    reads: &[EncodedRead],
    contigs: &Contigs,
) -> Option<Vec<CriticalRegion>> {
    if !contigs.is_repeat(from) {
        return None;
    }
    debug!("Contig:{} comes from repetitive contig!", from);
    let mut counts: Vec<u64> = vec![0; contigs.get_last_unit(from as u16)? as usize + 1];
    for read in reads {
        let mut prev_unit: Option<&Encode> = None;
        for unit in read.seq().iter().filter_map(|e| match e {
            ChunkedUnit::En(e) => Some(e),
            ChunkedUnit::Gap(_g) => None,
        }) {
            if let Some(prev_unit) = prev_unit {
                if prev_unit.contig != unit.contig && unit.contig == from {
                    // Entering the repetitive contig. Increment the landing point.
                    counts[unit.unit as usize] |= 1 << prev_unit.contig;
                } else if prev_unit.contig == from && unit.contig != from {
                    // Exiting the repetitive contig. Increment the leaving point.
                    counts[prev_unit.unit as usize] |= 1 << unit.contig;
                }
            }
            prev_unit = Some(unit);
        }
    }
    // Here, we ignore first 25 units, because it would be the 'canonical' connections.
    // TODO: Tune this.
    let len = counts.len();
    for i in 0..25.min(len) {
        counts[i] = 0;
        counts[len - 1 - i] = 0;
    }
    let counts: Vec<usize> = counts
        .into_iter()
        .map(|e| e.count_ones() as usize)
        .collect();
    debug!("Repetitive Summary");
    for (idx, count) in counts.iter().enumerate().filter(|&(_, &l)| l > 0) {
        debug!("{}\t{}", idx, count);
    }
    let region = calculate_average_more_than(&counts, 1);
    debug!("{:?}", region);
    let region = merge_overlap_and_remove_contained(region);
    debug!("{:?}", region);
    let result = region
        .into_iter()
        .map(|(start, end)| tighten_up(start, end, &counts))
        .map(|(start, end)| CriticalRegion::RJ(RepeatJunction::new(from, start as u16, end as u16)))
        .collect();
    Some(result)
}

fn calculate_average_more_than(input: &[usize], average: i32) -> Vec<(usize, usize)> {
    // Calculate maximal region with average more than READ_NUM.
    // 1. Extract average.
    let input: Vec<i32> = input.into_iter().map(|&e| e as i32 - average).collect();
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
    for i in 0..cumsum.len() {
        assert!(rmin[i] <= lmax[i]);
    }
    while i <= input.len() && j <= input.len() {
        if rmin[i] <= lmax[j] {
            j += 1;
        } else {
            // Output j and i.
            if 0 < j && rmin[i] <= lmax[j - 1] && i != j - 1 {
                result.push((i, j - 1));
            }
            while lmax[j] < rmin[i] {
                i += 1;
            }
        }
    }
    if 0 < j && rmin[i] <= lmax[j - 1] && i != j - 1 {
        result.push((i, j - 1));
    }
    result
}

fn tighten_up(start: usize, end: usize, xs: &[usize]) -> (usize, usize) {
    // Tighten up the region xs[start..end)
    // into xs[start+x..end-y), x and y are
    // the largest number that can increase the average score continuously.
    let mut average = xs[start..end].iter().sum::<usize>() as f64 / (end - start) as f64;
    // debug!("{}-{},average {}", start, end, average);
    let (mut start, mut end) = (start, end);
    while start + 1 < end {
        // if average would increase (ave < (ave * len-x)/(len-1)), decrease end by one.
        // the condition is the same as ...
        if (xs[end - 1] as f64) < average {
            end -= 1;
            average = xs[start..end].iter().sum::<usize>() as f64 / (end - start) as f64;
        } else {
            break;
        }
    }
    while start + 1 < end {
        if (xs[start] as f64) < average {
            start += 1;
            average = xs[start..end].iter().sum::<usize>() as f64 / (end - start) as f64;
        } else {
            break;
        }
    }
    assert!(start <= end);
    let _average = xs[start..end].iter().sum::<usize>() as f64 / (end - start) as f64;
    // debug!("Fin. {}-{},average {}", start, end, average);

    (start, end)
}

fn merge_overlap_and_remove_contained(mut regions: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if regions.is_empty() {
        return vec![];
    }
    regions.sort();
    let mut end = 0;
    let contained_free: Vec<_> = regions
        .into_iter()
        .filter(|&(_, e)| {
            if end < e {
                end = e;
                true
            } else {
                false
            }
        })
        .collect();
    let (mut start, mut end) = contained_free[0];
    let mut result = vec![];
    for &(s, e) in &contained_free[1..] {
        assert!(end < e);
        if end + 4 < s {
            // No overlap.
            result.push((start, end));
            start = s;
            end = e;
        } else {
            // Overlap.
            end = e;
        }
    }
    result.push((start, end));
    result
}

#[cfg(test)]
pub mod tests {
    use super::*;
    #[test]
    fn ave_max_test() {
        let input = vec![0, 0, 1, 1, 1, 0, 0, 0];
        let res = calculate_average_more_than(&input, 1);
        assert_eq!(res, vec![(2, 5)]);
    }
    #[test]
    fn ave_max_test2() {
        let input = vec![1, 1, 2, 1, 0, 0, 0, 0, 0];
        let res = calculate_average_more_than(&input, 1);
        assert_eq!(res, vec![(0, 5)]);
        let input = vec![1, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 1];
        let res = calculate_average_more_than(&input, 1);
        assert_eq!(res, vec![(0, 5), (9, 12)]);
    }
    #[test]
    fn tighen_up() {
        let input = vec![0, 0, 0, 0, 2, 1, 1, 1, 0];
        let (s, e) = tighten_up(3, 8, &input);
        assert_eq!((s, e), (4, 8));
    }
    #[test]
    fn get_max_min_unit_test() {
        let gen = |x, y| ChunkedUnit::En(Encode::sketch(x, y, true));
        let gen_gap = |x, y| ChunkedUnit::Gap(GapUnit::new(&vec![], Some((x, y))));
        let read = EncodedRead::from("test".to_string(), vec![gen(1, 2)]);
        assert_eq!(get_max_min_unit(&read, 0), (std::i32::MAX, std::i32::MIN));
        assert_eq!(get_max_min_unit(&read, 1), (2, 2));
        let read = EncodedRead::from(
            "test".to_string(),
            vec![gen(1, 2), gen(1, 3), gen(1, 4), gen(1, 1)],
        );
        assert_eq!(get_max_min_unit(&read, 1), (1, 4));
        let read = vec![
            gen(1, 2),
            gen(1, 3),
            gen(0, 121),
            gen(0, 89),
            gen(0, 12),
            gen(1, 0),
            gen(1, 21),
        ];
        let read = EncodedRead::from("test".to_string(), read);
        assert_eq!(get_max_min_unit(&read, 1), (0, 21));
        assert_eq!(get_max_min_unit(&read, 0), (12, 121));
        let read = vec![
            gen_gap(1, 10),
            gen(1, 10),
            gen(1, 11),
            gen(1, 12),
            gen(1, 13),
            gen_gap(10, 10),
            gen(1, 12),
            gen(1, 0),
        ];
        let read = EncodedRead::from("test".to_string(), read);
        assert_eq!(get_max_min_unit(&read, 1), (0, 13));
        let read = vec![gen_gap(10, 10), gen(10, 2), gen_gap(10, 11)];
        let read = EncodedRead::from("test".to_string(), read);
        assert_eq!(get_max_min_unit(&read, 10), (2, 2));
        assert_eq!(get_max_min_unit(&read, 0), (std::i32::MAX, std::i32::MIN));
    }
}
