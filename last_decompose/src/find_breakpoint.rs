use super::ERead;
use last_tiling::repeat::RepeatPairs;
use last_tiling::Contigs;
use serde::{Deserialize, Serialize};
use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
/// The threshould used to determine whether two regions would be regarded as pairs.
const READ_NUM: usize = 10;
/// How many reads are needed to construct a separate class.
/// Note that it should be more than 1, since there is a
/// danger of chimeric reads.
const CHECK_THR: u16 = 10;
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CriticalRegion {
    CP(ContigPair),
    CR(ConfluentRegion),
}

impl Display for CriticalRegion {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        match self {
            CriticalRegion::CP(ref cp) => write!(f, "{}", cp),
            CriticalRegion::CR(ref cr) => write!(f, "{}", cr),
        }
    }
}

// A trait to classify a given read.
// (read class, read index, read itself)
// type ReadWithClassAndIndex<'a> = (usize, usize, &'a EncodedRead);
pub trait ReadClassify {
    /// Whether the read spans `self`. In other words,
    /// If r is not agree with this critical region, return true.
    fn is_spanned_by(&self, r: &ERead) -> bool;
    /// Wether `self` contains the read.
    fn contains(&self, r: &ERead) -> bool;
    /// If r is along with this critical region, return true.
    fn along_with(&self, r: &ERead) -> bool;
}

impl ReadClassify for CriticalRegion {
    fn is_spanned_by(&self, r: &ERead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.is_spanned_by(r),
            CriticalRegion::CR(ref cr) => cr.is_spanned_by(r),
        }
    }
    fn contains(&self, r: &ERead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.contains(r),
            CriticalRegion::CR(ref cr) => cr.contains(r),
        }
    }
    fn along_with(&self, r: &ERead) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.along_with(r),
            CriticalRegion::CR(ref cr) => cr.along_with(r),
        }
    }
}

/// The position at contigs.
#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq, Hash)]
pub struct Position {
    contig: u16,
    start_unit: u16,
    end_unit: u16,
    direction: Direction,
}

impl Ord for Position {
    fn cmp(&self, other: &Self) -> Ordering {
        use Ordering::*;
        if self.contig < other.contig {
            return Less;
        } else if self.contig > other.contig {
            return Greater;
        }
        match (
            self.start_unit.cmp(&other.start_unit),
            self.end_unit.cmp(&other.end_unit),
        ) {
            (Less, _) | (Equal, Less) => return Less,
            (Greater, _) | (Equal, Greater) => return Greater,
            (Equal, Equal) => {}
        }
        self.direction.cmp(&other.direction)
    }
}

impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub enum Direction {
    UpStream,
    DownStream,
}

impl Position {
    fn new(contig: u16, s: u16, e: u16, direction: Direction) -> Self {
        let (start_unit, end_unit) = (s, e);
        Self {
            contig,
            start_unit,
            end_unit,
            direction,
        }
    }
    pub fn range(&self) -> (i32, i32) {
        (self.start_unit as i32, self.end_unit as i32)
    }
    pub fn contig(&self) -> u16 {
        self.contig
    }
    fn format(&self, name: &str) -> String {
        let (header, footer) = match self.direction {
            Direction::DownStream => ("     ", "<--->"),
            Direction::UpStream => ("<--->", "     "),
        };
        let s = format!("{:5}", self.start_unit);
        let t = format!("{:5}", self.end_unit);
        format!("{}{}:{}:{}{}", header, s, name, t, footer)
    }
    fn is_spanned_by(&self, r: &ERead) -> bool {
        let s = self.start_unit;
        let t = self.end_unit;
        let (s_thr, t_thr) = (
            s.max(CHECK_THR as u16) - CHECK_THR as u16,
            t + CHECK_THR as u16,
        );
        r.does_touch(self.contig, s_thr, s) && r.does_touch(self.contig, t, t_thr)
    }
    fn along_with(&self, r: &ERead) -> bool {
        let s = self.start_unit;
        let t = self.end_unit;
        let (s_thr, t_thr) = (
            s.max(CHECK_THR as u16) - CHECK_THR as u16,
            t + CHECK_THR as u16,
        );
        let c = self.contig;
        match self.direction {
            Direction::DownStream => !r.does_touch(c, s_thr, s) && r.does_touch(c, t, t_thr),
            Direction::UpStream => r.does_touch(c, s_thr, s) && !r.does_touch(c, t, t_thr),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq, Hash)]
pub struct ContigPair {
    contig1: Position,
    contig2: Position,
}

impl Display for ContigPair {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        writeln!(f, "{}", self.contig1.format("Recom"))?;
        write!(f, "{}", self.contig2.format("Recom"))
    }
}

impl Ord for ContigPair {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.contig1.cmp(&other.contig1) {
            Ordering::Equal => self.contig1.cmp(&other.contig2),
            x => x,
        }
    }
}
impl PartialOrd for ContigPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Return the maximum/minumum index of encoded read's unit with contig of `contig`.
/// Note that if the start unit is 0-th unit of `contig`, followed by large gap>UNIT_SIZE
/// It assigns -1.
/// And if the reads ended at i-th unit of `contig`, followed by large gap>UNIT_SIZE,
/// it assigns i+1.
/// If there's no unit in contig, return (std::i32::MAX, std::i32::MIN).
pub fn get_max_min_unit(r: &ERead, contig: u16) -> (i32, i32) {
    let filtered = r
        .seq()
        .iter()
        .filter(|c| c.contig() as u16 == contig)
        .map(|c| c.unit());
    let min = filtered.clone().min();
    let max = filtered.clone().max();
    match (min, max) {
        (Some(min), Some(max)) => (min as i32, max as i32),
        (None, None) => (std::i32::MAX, std::i32::MIN),
        _ => unreachable!(),
    }
}

impl ReadClassify for ContigPair {
    // The broad repeat is like the below.
    // -----A||||||B----
    // -----C||||||D----
    // Then, the read `r` spans when A->B or C->D.
    // Note that if r is not the same contig as self, or do not overlap at all,
    // it return false.
    fn is_spanned_by(&self, r: &ERead) -> bool {
        // Check if r spans A->B.
        let span_ab = self.contig1.is_spanned_by(r);
        // Check if r spans C->D.
        let span_cd = self.contig2.is_spanned_by(r);
        span_ab || span_cd
    }
    // Return whether self contains the read.
    // It is enough to check s < start && e < end holds
    // for contig 1 or contig 2.
    // Remember that, if there's no unit,
    // get_max_min_unit(,) would return (MAX,MIN), where
    // the condition above never holds.
    fn contains(&self, r: &ERead) -> bool {
        let (c1_s, c1_e) = self.contig1.range();
        let (c1_min, c1_max) = get_max_min_unit(r, self.contig1.contig);
        let contained_in_c1 = c1_s < c1_min && c1_max < c1_e;
        let (c2_s, c2_e) = self.contig1.range();
        let (c2_min, c2_max) = get_max_min_unit(r, self.contig2.contig);
        let contained_in_c2 = c2_s < c2_min && c2_max < c2_e;
        contained_in_c1 || contained_in_c2
    }
    // Check wether this read is along with self.
    // Note that this method is direction-sensitive.
    fn along_with(&self, r: &ERead) -> bool {
        self.contig1.along_with(r) && self.contig2.along_with(r)
    }
}

impl ContigPair {
    pub fn new(contig1: Position, contig2: Position) -> Self {
        Self { contig1, contig2 }
    }
    pub fn contig1(&self) -> &Position {
        &self.contig1
    }
    pub fn contig2(&self) -> &Position {
        &self.contig2
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq, Hash)]
pub struct ConfluentRegion {
    // Contains start and end position.
    pos: Position,
}

impl Display for ConfluentRegion {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        write!(f, "{}", self.pos.format("Entry"))
    }
}

impl Ord for ConfluentRegion {
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos.cmp(&other.pos)
    }
}
impl PartialOrd for ConfluentRegion {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ConfluentRegion {
    fn new(pos: Position) -> Self {
        Self { pos }
    }
    pub fn contig(&self) -> &Position {
        &self.pos
    }
}

impl ReadClassify for ConfluentRegion {
    // Check if `r` spans `self`. In other words,
    // It returns true if `r` do not agree with self.
    // Even thoguth this function return false,
    // It does not garantee that `r` is along with this junction.
    fn is_spanned_by(&self, r: &ERead) -> bool {
        self.pos.is_spanned_by(r)
    }
    fn contains(&self, r: &ERead) -> bool {
        let (min, max) = get_max_min_unit(r, self.pos.contig);
        let (start, end) = self.pos.range();
        start < min && max < end
    }
    fn along_with(&self, r: &ERead) -> bool {
        self.pos.along_with(r)
    }
}

/// Return critical regions.
pub fn critical_regions(
    reads: &[ERead],
    contigs: &Contigs,
    repeats: &[RepeatPairs],
) -> Vec<CriticalRegion> {
    let num_of_contig = contigs.get_num_of_contigs() as u16;
    let mut regions = vec![];
    debug!("There are {} contigs", num_of_contig);
    debug!("There are {} reads", reads.len());
    for c in contigs.names().iter() {
        debug!("->{}({}len)", c, contigs.get(c).unwrap().len());
    }
    for from in 0..num_of_contig {
        if let Some(res) = critical_region_within(from, reads, contigs, repeats) {
            // for c in &res {
            //     debug!("{:?}", c);
            // }
            regions.extend(res);
        }
        debug!("Inner end");
        debug!("Confluent region");
        if let Some(res) = critical_region_confluent(from, reads, contigs) {
            // for c in &res {
            //     debug!("{:?}", c);
            // }
            regions.extend(res);
        }
    }
    regions
}

// Return all critical region within a given contig.
// Note that the confluent region would not be captured.
fn critical_region_within(
    contig: u16,
    reads: &[ERead],
    contigs: &Contigs,
    repeats: &[RepeatPairs],
) -> Option<Vec<CriticalRegion>> {
    debug!("{}-th contig, self critical region.", contig);
    let last_unit = contigs.get_last_unit(contig)? as usize;
    let mut inner_count: Vec<Vec<&ERead>> = vec![vec![]; last_unit + 1];
    const ACCEPTED_GAP: u16 = 10;
    for read in reads.iter().filter(|r| r.has(contig)) {
        for w in read.seq().windows(2) {
            if w[0].contig == contig && w[1].contig == contig {
                let gapsize = w[1].unit.max(w[0].unit) - w[1].unit.min(w[0].unit);
                let is_circled = last_unit as u16 - w[1].unit.max(w[0].unit) < ACCEPTED_GAP
                    && w[1].unit.min(w[0].unit) < ACCEPTED_GAP;
                if gapsize > ACCEPTED_GAP && !is_circled {
                    inner_count[w[0].unit as usize].push(read);
                    inner_count[w[1].unit as usize].push(read);
                }
            }
        }
    }
    // debug!("Profiled!({}len)", inner_count.len());
    // debug!("Position\tCount");
    // for (idx, count) in inner_count
    //     .iter()
    //     .map(|e| e.len())
    //     .enumerate()
    //     .filter(|&(_, len)| len != 0)
    // {
    //     debug!("{}\t{}", idx, count);
    // }
    let inner_region = {
        let inner_counts = inner_count.iter().map(|e| e.len()).collect::<Vec<_>>();
        let inner_region: Vec<_> = peak_detection(&inner_counts, READ_NUM)
            .into_iter()
            .chain(repeats.iter().flat_map(|repeats| {
                (repeats
                    .inner()
                    .iter()
                    .map(|r| (r.start_in_unit() as usize, r.end_in_unit() as usize)))
            }))
            .collect();
        merge_overlap_and_remove_contained(inner_region)
    };
    // debug!("Aggregated!");
    // debug!("Start\tEnd");
    // for &(s, t) in inner_region.iter() {
    //     debug!("{}\t{}", s, t);
    // }
    let mut result = vec![];
    for i in 0..inner_region.len() {
        for j in (i + 1)..inner_region.len() {
            let &(s, t) = &inner_region[i];
            let &(x, y) = &inner_region[j];
            let is_ordered = 0 < t && t <= x && x <= y;
            if !is_ordered {
                // (s,t) contains (x,y) or vise versa.
                continue;
            } else if x - t < 100 {
                // Too near.
                continue;
            }
            // Let's move on to case exhaution.
            // ----s||||t-----
            // ----x||||y-----
            // Case1: --->s(jumps)y---> : (s,t):UpStream, (x,y):Downstream
            // Case2: --->x(jumps)t---> : (s,t):DownStream, (x,y):UpStream
            // Case3: --->s(jumps)x---> : (s,t):UpStream, (x,y):UpStream,
            // Case4: --->t(jumps)y---> : (s,t):DownStream, (x,y):DownStream
            // To this end, first create each gadget.
            // Note that all the read innner_count should be mapped to the
            // 'contig' contig. However, to treat chimera reads correctly,
            // we double check it.
            // debug!("S-T:{}-{}\tX-Y:{}-{}", s, t, x, y);
            let (s_thr, t_thr) = (s.max(CHECK_THR) - CHECK_THR, t + CHECK_THR);
            let st_upstream: HashSet<_> = reads
                .iter()
                .filter(|read| {
                    read.does_touch(contig, s_thr, s) && !read.does_touch(contig, t, t_thr)
                })
                .collect();
            let st_downstream: HashSet<_> = reads
                .iter()
                .filter(|read| {
                    !read.does_touch(contig, s_thr, s) && read.does_touch(contig, t, t_thr)
                })
                .collect();
            let (x_thr, y_thr) = (x.max(CHECK_THR) - CHECK_THR, y + CHECK_THR);
            let xy_upstream: HashSet<_> = reads
                .iter()
                .filter(|read| {
                    read.does_touch(contig, x_thr, x) && !read.does_touch(contig, y, y_thr)
                })
                .collect();
            let xy_downstream: HashSet<_> = reads
                .iter()
                .filter(|read| {
                    !read.does_touch(contig, x, x_thr) && read.does_touch(contig, y, y_thr)
                })
                .collect();
            assert!(st_upstream.intersection(&st_downstream).count() == 0);
            assert!(xy_upstream.intersection(&xy_downstream).count() == 0);
            // Case1.
            if st_upstream.intersection(&xy_downstream).count() >= READ_NUM {
                debug!(
                    "Up-Down:{}",
                    st_upstream.intersection(&xy_downstream).count()
                );
                let c1 = Position::new(contig, s, t, Direction::UpStream);
                let c2 = Position::new(contig, x, y, Direction::DownStream);
                result.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case2
            if st_downstream.intersection(&xy_upstream).count() >= READ_NUM {
                debug!(
                    "Down-Up:{}",
                    st_downstream.intersection(&xy_upstream).count()
                );
                let c1 = Position::new(contig, s, t, Direction::DownStream);
                let c2 = Position::new(contig, x, y, Direction::UpStream);
                result.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case3
            if st_upstream.intersection(&xy_upstream).count() >= READ_NUM {
                let c1 = Position::new(contig, s, t, Direction::UpStream);
                let c2 = Position::new(contig, x, y, Direction::UpStream);
                result.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case4
            if st_downstream.intersection(&xy_downstream).count() >= READ_NUM {
                let c1 = Position::new(contig, s, t, Direction::DownStream);
                let c2 = Position::new(contig, x, y, Direction::DownStream);
                result.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
        }
    }
    Some(result)
}

#[allow(dead_code)]
// Currently we do not need this function. TODO.
fn critical_regions_between(
    from: u16,
    to: u16,
    reads: &[ERead],
    contigs: &Contigs,
) -> Option<Vec<CriticalRegion>> {
    // Enumerate critical regions from `from` contig to `to` contig.
    warn!("Invoke critical_regions_between. Currently not implemented correctly. Be careful.");
    let from_last_unit = contigs.get_last_unit(from)?;
    let to_last_unit = contigs.get_last_unit(to)?;
    let mut from_count: Vec<Vec<&ERead>> = vec![vec![]; from_last_unit as usize + 1];
    let mut to_count: Vec<Vec<&ERead>> = vec![vec![]; to_last_unit as usize + 1];
    for read in reads
        .iter()
        .filter(|r| r.has(from) && r.has(to) && r.len() > 0)
    {
        let mut prev = &read.seq()[0];
        for unit in read.seq().iter().skip(1) {
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
    // debug!("btw {}-{}", from, to);
    let from_region: Vec<_> = {
        let raw_count: Vec<_> = from_count.iter().map(|e| e.len()).collect();
        // debug!("Row count of {}", from);
        // for (idx, count) in raw_count
        //     .iter()
        //     .enumerate()
        //     .filter(|&(_, &count)| count > 0)
        // {
        //     debug!("{}\t{}", idx, count);
        // }
        let res = calculate_average_more_than(&raw_count, READ_NUM as i32);
        // for &(s, e) in &res {
        //     debug!("FROM:[{},{})", s, e);
        // }
        let res = merge_overlap_and_remove_contained(res);
        res.into_iter()
            //.map(|(s, e)| tighten_up(s, e, &raw_count))
            .map(|(s, e)| (s as u16, e as u16))
            .collect()
    };
    let to_region: Vec<_> = {
        let raw_count: Vec<_> = to_count.iter().map(|e| e.len()).collect();
        // debug!("Row count of {}", to);
        // for (idx, count) in raw_count
        //     .iter()
        //     .enumerate()
        //     .filter(|&(_, &count)| count > 0)
        // {
        //     debug!("{}\t{}", idx, count);
        // }
        let res = calculate_average_more_than(&raw_count, READ_NUM as i32);
        // for &(s, e) in &res {
        //     debug!("TO:[{},{})", s, e);
        // }
        let res = merge_overlap_and_remove_contained(res);
        res.into_iter()
            // .map(|(s, e)| tighten_up(s, e, &raw_count))
            .map(|(s, e)| (s as u16, e as u16))
            .collect()
    };
    // debug!("{:?}", from_region);
    // debug!("{:?}", to_region);
    let mut regions = vec![];
    for &(s, t) in &from_region {
        for &(x, y) in &to_region {
            // Let's move on to case exhaution.
            // From:----s||||t-----
            // To  :----x||||y-----
            // Case1: --->s(jumps)y---> : (s,t):UpStream, (x,y):Downstream
            // Case2: --->x(jumps)t---> : (s,t):DownStream, (x,y):UpStream
            // Case3: --->s(jumps)x---> : (s,t):UpStream, (x,y):UpStream,
            // Case4: --->t(jumps)y---> : (s,t):DownStream, (x,y):DownStream
            // To this end, first create each gadget.
            let st_upstream = from_count[s as usize..t as usize]
                .iter()
                .flat_map(|reads| {
                    reads.iter().filter(|read| {
                        let (min, max) = get_max_min_unit(read, from);
                        min <= s as i32 && s as i32 <= max
                    })
                })
                .copied();
            let st_downstream = from_count[s as usize..t as usize]
                .iter()
                .flat_map(|reads| {
                    reads.iter().filter(|read| {
                        let (min, max) = get_max_min_unit(read, from);
                        min <= t as i32 && t as i32 <= max
                    })
                })
                .copied();
            let xy_upstream = to_count[x as usize..y as usize]
                .iter()
                .flat_map(|reads| {
                    reads.iter().filter(|read| {
                        let (min, max) = get_max_min_unit(read, to);
                        min <= x as i32 && y as i32 <= max
                    })
                })
                .copied();
            let xy_downstream = to_count[x as usize..y as usize]
                .iter()
                .flat_map(|reads| {
                    reads.iter().filter(|read| {
                        let (min, max) = get_max_min_unit(read, to);
                        min <= y as i32 && y as i32 <= max
                    })
                })
                .copied();
            let thr = READ_NUM * 2;
            // Case1.
            let up_down: HashSet<&ERead> =
                st_upstream.clone().chain(xy_downstream.clone()).collect();
            if up_down.len() >= thr {
                let c1 = Position::new(from, s, t, Direction::UpStream);
                let c2 = Position::new(to, x, y, Direction::DownStream);
                regions.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case2
            let down_up: HashSet<&ERead> =
                st_downstream.clone().chain(xy_upstream.clone()).collect();
            if down_up.len() >= thr {
                let c1 = Position::new(from, s, t, Direction::DownStream);
                let c2 = Position::new(to, x, y, Direction::UpStream);
                regions.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case3
            let up_up: HashSet<&ERead> = st_upstream.chain(xy_upstream).collect();
            if up_up.len() >= thr {
                let c1 = Position::new(from, s, t, Direction::UpStream);
                let c2 = Position::new(to, x, y, Direction::UpStream);
                regions.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
            // Case4
            let down_down: HashSet<&ERead> = st_downstream.chain(xy_downstream).collect();
            if down_down.len() >= thr {
                let c1 = Position::new(from, s, t, Direction::DownStream);
                let c2 = Position::new(to, x, y, Direction::DownStream);
                regions.push(CriticalRegion::CP(ContigPair::new(c1, c2)));
            }
        }
    }
    Some(regions)
}

fn critical_region_confluent(
    contig: u16,
    reads: &[ERead],
    contigs: &Contigs,
) -> Option<Vec<CriticalRegion>> {
    // debug!("{}-th contig, self critical region.", contig);
    let last_unit = contigs.get_last_unit(contig)? as usize;
    let mut clip_count: Vec<Vec<&ERead>> = vec![vec![]; last_unit + 1];
    for read in reads.iter().filter(|r| r.has(contig) && r.len() > 2) {
        let first_chunk = &read.seq()[0];
        if first_chunk.contig == contig {
            clip_count[first_chunk.unit as usize].push(read);
        }
        let last_chunk = &read.seq().last().unwrap();
        if last_chunk.contig == contig && last_chunk.unit != last_unit as u16 {
            clip_count[last_chunk.unit as usize].push(read);
        }
    }
    // debug!("Profiled!({}len)", clip_count.len());
    // debug!("Position\tCount");
    // for (idx, count) in clip_count
    //     .iter()
    //     .map(|e| e.len())
    //     .enumerate()
    //     .filter(|&(_, len)| len != 0)
    // {
    //     debug!("{}\t{}", idx, count);
    // }
    let clip_region = {
        let clip_counts = clip_count.iter().map(|e| e.len()).collect::<Vec<_>>();
        let mean = reads.len() / last_unit * 2;
        let clip_region = peak_detection(&clip_counts, READ_NUM + 2 * mean);
        merge_overlap_and_remove_contained(clip_region)
            .into_iter()
            .map(|(s, e)| (s as u16, e as u16))
            .collect::<Vec<_>>()
    };
    // debug!("Aggregated!");
    // debug!("Start\tEnd");
    // for &(s, t) in clip_region.iter() {
    //     debug!("{}\t{}", s, t);
    // }
    let mut result = vec![];
    for (s, t) in clip_region {
        if s.max(t) - s.min(t) > 20 {
            // Too long to be conlfulent regions.
            continue;
        }
        // Let's move on to case exhaution.
        // ----s||||t-----
        // Case1: <---s(entry)y     : (s,t):UpStream
        // Case2:     s(entry)t---> : (s,t):DownStream
        let (s_thr, t_thr) = (s.max(CHECK_THR) - CHECK_THR, t + CHECK_THR);
        // Here, we maybe want to get circulaized genome region.
        // It should not happen in x-y, because we force 0 < t < x < y.
        let st_upstream: HashSet<&ERead> = reads
            .iter()
            .filter(|read| read.does_touch(contig, s_thr, s) && !read.does_touch(contig, t_thr, t))
            .collect();
        let st_downstream: HashSet<&ERead> = reads
            .iter()
            .filter(|read| !read.does_touch(contig, s_thr, s) && read.does_touch(contig, t_thr, t))
            .collect();
        let thr = READ_NUM + 2 * reads.len() / last_unit;
        // Case1.
        if st_upstream.len() >= thr {
            let c1 = Position::new(contig, s, t, Direction::UpStream);
            result.push(CriticalRegion::CR(ConfluentRegion::new(c1)));
        }
        // Case2
        if st_downstream.len() >= thr {
            let c1 = Position::new(contig, s, t, Direction::DownStream);
            result.push(CriticalRegion::CR(ConfluentRegion::new(c1)));
        }
    }
    Some(result)
}

fn calculate_average_more_than(input: &[usize], average: i32) -> Vec<(usize, usize)> {
    // Calculate maximal region with average more than READ_NUM.
    // 1. Extract average.
    let input: Vec<i32> = input.iter().map(|&e| e as i32 - average).collect();
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

fn merge_overlap_and_remove_contained(mut regions: Vec<(usize, usize)>) -> Vec<(u16, u16)> {
    if regions.is_empty() {
        return vec![];
    }
    regions.sort();
    let (mut start, mut end) = regions[0];
    let mut result = vec![];
    for &(s, e) in &regions[1..] {
        if end + 1 < s {
            // No overlap.
            result.push((start, end));
            start = s;
            end = e;
        } else {
            // Overlap.
            end = e.max(end);
        }
    }
    result.push((start, end));
    result
        .into_iter()
        .map(|(s, t)| (s as u16, t as u16))
        .collect()
}

fn peak_detection(counts: &[usize], threshold: usize) -> Vec<(usize, usize)> {
    let (res, _, _) = counts.iter().enumerate().fold(
        (vec![], false, 0),
        |(mut res, is_in_peak, start), (idx, &count)| match (is_in_peak, count < threshold) {
            (true, true) => {
                res.push((start, idx));
                (res, false, idx)
            }
            (true, false) => (res, true, start),
            (false, true) => (res, false, idx),
            (false, false) => (res, true, idx),
        },
    );
    res
}

#[allow(dead_code)]
fn merge_circular_genome(mut regions: Vec<(u16, u16)>, len: usize) -> Vec<(u16, u16)> {
    let len = len as u16;
    if regions.len() <= 1 {
        regions
    } else {
        let first: (u16, u16) = *regions.first().unwrap();
        let last: (u16, u16) = *regions.last().unwrap();
        if first.0 < 10 && last.1 - len < 10 {
            regions.remove(0);
            regions.pop();
            regions.insert(0, (last.0, first.1));
        }
        regions
    }
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
