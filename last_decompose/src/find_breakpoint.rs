#![allow(dead_code)]
use super::ERead;
use last_tiling::repeat::RepeatPairs;
use last_tiling::{Contigs, LastTAB};
use serde::{Deserialize, Serialize};
use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
const MERGE_THR: usize = 500;
const COVERAGE_THR: usize = 25;
const WINDOW: usize = 3;
use std::collections::HashMap;

// The threshould used to determine whether two regions would be regarded as pairs.
const fn factor(x: usize) -> usize {
    x * 3 / 2
}
/// How many reads are needed to construct a separate class.
/// Note that it should be more than 1, since there is a
/// danger of chimeric reads.
const CHECK_THR: u16 = 2;
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
    longest: u16,
}

impl Ord for Position {
    fn cmp(&self, other: &Self) -> Ordering {
        use Ordering::*;
        match self.contig.cmp(&other.contig) {
            Equal => {}
            x => return x,
        };
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
    fn new(contig: u16, s: u16, e: u16, direction: Direction, max: u16) -> Self {
        let (start_unit, end_unit) = (s, e);
        Self {
            contig,
            start_unit,
            end_unit,
            direction,
            longest: max,
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
        let (s_thr, t_thr) = (s.max(CHECK_THR) - CHECK_THR, t + CHECK_THR);
        if t_thr < self.longest {
            r.does_touch(self.contig, s_thr, s) && r.does_touch(self.contig, t, t_thr)
        } else {
            r.does_touch(self.contig, s_thr, s) && r.does_touch(self.contig, 0, CHECK_THR)
        }
    }
    fn along_with(&self, r: &ERead) -> bool {
        let s = self.start_unit;
        let t = self.end_unit;
        let (s_thr, t_thr) = (s.max(CHECK_THR) - CHECK_THR, t + CHECK_THR);
        let c = self.contig;
        match self.direction {
            Direction::UpStream => r.does_touch(c, s_thr, s) && !r.does_touch(c, t, t_thr),
            Direction::DownStream => !r.does_touch(c, s_thr, s) && r.does_touch(c, t, t_thr),
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
    _repeats: &[RepeatPairs],
    alignments: &[LastTAB],
) -> Vec<CriticalRegion> {
    let num_of_contig = contigs.get_num_of_contigs() as u16;
    debug!("There are {} contigs", num_of_contig);
    debug!("There are {} reads", reads.len());
    for c in contigs.names().iter() {
        debug!("->{}({}len)", c, contigs.get(c).unwrap().len());
    }
    let contig_pairs = contigpair_position(reads, contigs);
    let confluent_regions = confluent_position(alignments, contigs, last_tiling::UNIT_SIZE);
    let contig_pairs = contig_pairs.into_iter().map(|cp| CriticalRegion::CP(cp));
    let confluent_regions = confluent_regions
        .into_iter()
        .map(|cr| CriticalRegion::CR(cr));
    contig_pairs.chain(confluent_regions).collect()
}

pub fn contigpair_position(reads: &[ERead], contigs: &Contigs) -> Vec<ContigPair> {
    let mut from_up_counts: Vec<Vec<Vec<(u16, usize, bool)>>> = contigs
        .get_last_units()
        .iter()
        .map(|&len| vec![vec![]; len as usize + 1])
        .collect();
    let max_len: Vec<_> = contigs
        .get_last_units()
        .into_iter()
        .map(|e| e as usize)
        .collect();
    let mut to_down_counts = from_up_counts.clone();
    for read in reads.iter() {
        let jumps = read.seq().windows(4).filter(is_jumping);
        for w in jumps {
            let (start, end) = (w[1].unit as usize, w[2].unit as usize);
            let (c1, c2) = (w[1].contig as usize, w[2].contig as usize);
            // Check if they are not circular.
            if c1 == c2 && (start.min(end) < 5 && max_len[c1 as usize] - 5 < end.max(start)) {
                continue;
            }
            // If w[0] < w[1], it is to upstream.
            let former_from_up = is_from_up(&w[..2]);
            // If w[2] < w[3], it is to downstream.
            let latter_from_up = !is_from_up(&w[2..]);
            match (start < end, former_from_up, latter_from_up) {
                (true, true, x) => from_up_counts[c1][start].push((c2 as u16, end, x)),
                (true, false, x) => to_down_counts[c1][start].push((c2 as u16, end, x)),
                (false, x, true) => from_up_counts[c2][end].push((c1 as u16, start, x)),
                (false, x, false) => to_down_counts[c2][end].push((c1 as u16, start, x)),
            }
        }
    }
    for contig in from_up_counts.iter() {
        for (idx, jumps) in contig.iter().enumerate() {
            if jumps.len() > 2 {
                debug!("U\t{}\t{}", idx, jumps.len());
            }
        }
    }
    for contig in to_down_counts.iter() {
        for (idx, jumps) in contig.iter().enumerate() {
            if jumps.len() > 2 {
                debug!("D\t{}\t{}", idx, jumps.len());
            }
        }
    }
    let mut contigpairs = vec![];
    //debug!("start call peak");
    let pairs = peak_call_contigpair(from_up_counts, true, contigs);
    contigpairs.extend(pairs);
    //debug!("start call peak");
    let pairs = peak_call_contigpair(to_down_counts, false, contigs);
    contigpairs.extend(pairs);
    contigpairs
}

fn peak_call_contigpair(
    mut counts: Vec<Vec<Vec<(u16, usize, bool)>>>,
    from_upstream: bool,
    contigs: &Contigs,
) -> Vec<ContigPair> {
    let mut contigpairs = vec![];
    for (contig, jumps) in counts.iter_mut().enumerate() {
        let from_max = contigs.get_last_unit(contig as u16).unwrap();
        let mut start_idx = 0;
        while start_idx < jumps.len() {
            // Search a heviest edge from this contig.
            let (from, to_c, to, direction) = match search_jump_start(jumps, start_idx) {
                Some(res) => res,
                None => break,
            };
            //debug!("{}:{}->{}:{}({})", contig, from, to_c, to, direction);
            // Determine the range of the start position ad the range of
            // the end position.
            let (from_start, from_end, to_start, to_end) =
                search_jump_area(jumps, from, to_c, to, direction);
            // debug!("{}-{} -> {}-{}", from_start, from_end, to_start, to_end);
            // Collect the edges from (start,end) to (start_dst, end_dst)
            let counts = collect_edges_in(
                jumps, from_start, from_end, to_c, to_start, to_end, direction,
            );
            // debug!("having {} reads", counts);
            // and remove them from counts.
            if counts > COVERAGE_THR {
                remove_edges(
                    jumps, from_start, from_end, to_c, to_start, to_end, direction,
                );
                use Direction::*;
                let from_direction = if from_upstream { UpStream } else { DownStream };
                let to_direction = if direction { UpStream } else { DownStream };
                let to_max = contigs.get_last_unit(to_c as u16).unwrap();
                let (contig, from_start, from_end) =
                    (contig as u16, from_start as u16, from_end as u16);
                let pos1 = Position::new(contig, from_start, from_end, from_direction, from_max);
                let (to_start, to_end) = (to_start as u16, to_end as u16);
                let pos2 = Position::new(to_c, to_start, to_end, to_direction, to_max);
                contigpairs.push(ContigPair::new(pos1, pos2));
            }
            start_idx = from_end;
        }
    }
    contigpairs
}

const JUMP_THR: usize = 5;
// Return the start position of the edge, the contig of the destination, the position of the destimation, and whether the edge is 'to_downstream'.
fn search_jump_start(
    jumps: &[Vec<(u16, usize, bool)>],
    s: usize,
) -> Option<(usize, u16, usize, bool)> {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for (position, edges) in jumps.iter().enumerate().skip(s) {
        for key in edges.iter() {
            *counts.entry(key.clone()).or_default() += 1;
        }
        if let Some((&(x, y, z), _count)) =
            counts.iter().filter(|&(_, &count)| count > JUMP_THR).nth(0)
        {
            return Some((position, x, y, z));
        }
        counts.clear();
    }
    None
}

const OFFSET: usize = 5;
// Return the area of the junctions.
fn search_jump_area(
    jumps: &[Vec<(u16, usize, bool)>],
    from: usize,
    to_c: u16,
    to: usize,
    direction: bool,
) -> (usize, usize, usize, usize) {
    // Count the number of edges from from_start position to to_c:[to_start..to_end).
    let count = |pos: usize, start: usize, end: usize| {
        let start = start.max(OFFSET) - OFFSET;
        let end = end + OFFSET;
        jumps[pos]
            .iter()
            .filter(|&&(c, t, d)| c == to_c && start <= t && t < end && d == direction)
            .count()
    };
    // Get the minimum/maximum position from 'pos' to the target range.
    let get_range = |pos: usize, start: usize, end: usize| {
        let start = start.max(OFFSET) - OFFSET;
        let end = end + OFFSET;
        let edges = jumps[pos]
            .iter()
            .filter(|&&(c, t, d)| c == to_c && start <= t && t < end && d == direction);
        let min = edges
            .clone()
            .map(|&(_, t, _)| t)
            .min()
            .unwrap_or_else(|| panic!("{}", line!()));
        let max = edges
            .clone()
            .map(|&(_, t, _)| t)
            .max()
            .unwrap_or_else(|| panic!("{}", line!()));
        (min, max)
    };
    let (mut from_start, mut from_end) = (from, from + 1);
    let (mut to_start, mut to_end) = (to, to + 1);
    const STEP_SIZE: usize = 5;
    while let Some(step) = (1..STEP_SIZE)
        .filter(|&n| n < from_start && count(from_start - n, to_start, to_end) > 0)
        .last()
    {
        from_start -= step;
        let (start, end) = get_range(from_start, to_start, to_end);
        to_start = to_start.min(start);
        to_end = to_end.max(end);
    }
    while let Some(step) = (1..STEP_SIZE)
        .filter(|&n| n + from_end <= jumps.len() && count(from_end + n - 1, to_start, to_end) > 0)
        .last()
    {
        from_end += step;
        let (start, end) = get_range(from_end - 1, to_start, to_end);
        to_start = to_start.min(start);
        to_end = to_end.max(end);
    }
    (from_start, from_end, to_start, to_end)
}

// Count the number of edges from [f_start, f_end) and ending in to_c:[to_start, to_end).
fn collect_edges_in(
    jumps: &[Vec<(u16, usize, bool)>],
    f_start: usize,
    f_end: usize,
    to_c: u16,
    to_start: usize,
    to_end: usize,
    direction: bool,
) -> usize {
    jumps[f_start..f_end]
        .iter()
        .map(|edges| {
            edges
                .iter()
                .filter(|&&(c, t, d)| c == to_c && d == direction && to_start <= t && t < to_end)
                .count()
        })
        .sum::<usize>()
}

fn remove_edges(
    jumps: &mut [Vec<(u16, usize, bool)>],
    f_start: usize,
    f_end: usize,
    to_c: u16,
    to_start: usize,
    to_end: usize,
    direction: bool,
) {
    for edges in jumps[f_start..f_end].iter_mut() {
        *edges = edges
            .iter()
            .filter(|&&(c, t, d)| !(c == to_c && d == direction && to_start <= t && t < to_end))
            .copied()
            .collect();
    }
}

fn is_jumping(w: &&[super::CUnit]) -> bool {
    const ACCEPTED_GAP: u16 = 10;
    fn diff(x: &super::CUnit, y: &super::CUnit) -> u16 {
        (x.unit).max(y.unit) - (x.unit).min(y.unit)
    }
    let is_first_part_cont = w[0].contig == w[1].contig && diff(&w[0], &w[1]) <= ACCEPTED_GAP;
    let is_second_part_cont = w[2].contig == w[3].contig && diff(&w[2], &w[3]) <= ACCEPTED_GAP;
    let is_cont = is_first_part_cont && is_second_part_cont;
    is_cont && (w[1].contig != w[2].contig || diff(&w[1], &w[2]) > ACCEPTED_GAP)
}

fn is_from_up(w: &[super::CUnit]) -> bool {
    w[0].unit < w[1].unit
}

// Return all critical region within a given contig.
// Note that the confluent region would not be captured.
fn critical_region_within(
    contig: u16,
    reads: &[ERead],
    contigs: &Contigs,
    repeats: &[RepeatPairs],
) -> Option<Vec<ContigPair>> {
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
    let ave = reads.len() / last_unit;
    let inner_region = {
        let inner_counts = inner_count.iter().map(|e| e.len()).collect::<Vec<_>>();
        let inner_region: Vec<_> = peak_detection(&inner_counts, factor(ave))
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
    let mut result = vec![];
    for i in 0..inner_region.len() {
        for j in (i + 1)..inner_region.len() {
            let &(s, t) = &inner_region[i];
            let &(x, y) = &inner_region[j];
            let is_ordered = 0 < t && t <= x && x <= y;
            if !is_ordered || x - t < 100 {
                // (s,t) contains (x,y) or vise versa.
                // Or, Too near.
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
            let width = (t - s).max(y - x) as usize;
            let thr = factor(width * ave);
            // Case1.
            let last_unit = last_unit as u16;
            if st_upstream.intersection(&xy_downstream).count() >= thr {
                let c1 = Position::new(contig, s, t, Direction::UpStream, last_unit);
                let c2 = Position::new(contig, x, y, Direction::DownStream, last_unit);
                result.push(ContigPair::new(c1, c2));
            }
            // Case2
            if st_downstream.intersection(&xy_upstream).count() >= thr {
                // debug!(
                //     "Down-Up:{}",
                //     st_downstream.intersection(&xy_upstream).count()
                // );
                let c1 = Position::new(contig, s, t, Direction::DownStream, last_unit);
                let c2 = Position::new(contig, x, y, Direction::UpStream, last_unit);
                result.push(ContigPair::new(c1, c2));
            }
            // Case3
            if st_upstream.intersection(&xy_upstream).count() >= thr {
                let c1 = Position::new(contig, s, t, Direction::UpStream, last_unit);
                let c2 = Position::new(contig, x, y, Direction::UpStream, last_unit);
                result.push(ContigPair::new(c1, c2));
            }
            // Case4
            if st_downstream.intersection(&xy_downstream).count() >= thr {
                let c1 = Position::new(contig, s, t, Direction::DownStream, last_unit);
                let c2 = Position::new(contig, x, y, Direction::DownStream, last_unit);
                result.push(ContigPair::new(c1, c2));
            }
        }
    }
    Some(result)
}

fn merge_overlap_and_remove_contained(mut regions: Vec<(usize, usize)>) -> Vec<(u16, u16)> {
    if regions.is_empty() {
        return vec![];
    }
    regions.sort();
    let (mut start, mut end) = regions[0];
    let mut result = vec![];
    for &(s, e) in &regions[1..] {
        if end + 10 < s {
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

/// Enumerate confluent region in the references.
pub fn confluent_position(
    alignments: &[LastTAB],
    contigs: &Contigs,
    unit_size: usize,
) -> Vec<ConfluentRegion> {
    let references: Vec<_> = {
        let mut reference: Vec<_> = alignments
            .iter()
            .map(|r| (r.seq1_name(), r.seq1_len()))
            .collect();
        reference.sort_by_key(|&e| e.1);
        reference.dedup();
        reference
    };
    let mut start_stop_count: HashMap<String, Vec<(usize, usize)>> = references
        .into_iter()
        .map(|(id, len)| (id.to_string(), vec![(0, 0); len + 1]))
        .collect();
    let mut reads: HashMap<_, Vec<_>> = HashMap::new();
    for tab in alignments {
        reads
            .entry(tab.seq2_name().to_string())
            .or_insert_with(Vec::new)
            .push(tab);
    }
    for (_id, tabs) in reads.into_iter() {
        if let Some((ref_name, position)) = get_first_alignment(&tabs) {
            if let Some(res) = start_stop_count.get_mut(ref_name) {
                res[position].0 += 1;
            }
        }
        if let Some((ref_name, position)) = get_last_alignment(&tabs) {
            if let Some(res) = start_stop_count.get_mut(ref_name) {
                res[position].1 += 1;
            }
        }
    }
    let mut result = vec![];
    for (id, counts) in start_stop_count {
        let downstream: Vec<_> = counts.iter().map(|x| x.0).collect();
        let chunks = peak_call(&downstream)
            .into_iter()
            .map(|(s, e)| (id.clone(), s, e, true));
        result.extend(chunks);
        let upstream: Vec<_> = counts.iter().map(|x| x.1).collect();
        let chunks = peak_call(&upstream)
            .into_iter()
            .map(|(s, e)| (id.clone(), s, e, false));
        result.extend(chunks);
    }
    result
        .iter()
        .map(|&(ref id, start, end, to_downstream)| {
            let id = contigs.get_id(id).unwrap();
            let (start, end) = (start / unit_size, end / unit_size);
            let (start, end) = (start as u16, end as u16);
            let max = contigs.get_last_unit(id).unwrap();
            use Direction::*;
            let di = if to_downstream { DownStream } else { UpStream };
            ConfluentRegion::new(Position::new(id, start, end, di, max))
        })
        .collect()
}

fn get_first_alignment<'a>(tabs: &[&'a LastTAB]) -> Option<(&'a str, usize)> {
    let tab = tabs.iter().min_by_key(|t| t.seq2_start_from_forward())?;
    let id = tab.seq1_name();
    let position = tab.seq1_start_from_forward();
    Some((id, position))
}

fn get_last_alignment<'a>(tabs: &[&'a LastTAB]) -> Option<(&'a str, usize)> {
    let tab = tabs.iter().min_by_key(|t| t.seq2_end_from_forward())?;
    let id = tab.seq1_name();
    let position = tab.seq1_end_from_forward();
    Some((id, position))
}

fn peak_call(counts: &[usize]) -> Vec<(usize, usize)> {
    let positions: Vec<_> = counts
        .windows(WINDOW)
        .enumerate()
        .filter(|(_, w)| w.iter().sum::<usize>() >= COVERAGE_THR)
        .map(|(idx, _)| idx)
        .collect();
    merge(positions, MERGE_THR)
}

fn merge(positions: Vec<usize>, thr: usize) -> Vec<(usize, usize)> {
    let mut result: Vec<(usize, usize)> = vec![];
    let mut start = 0;
    while start < positions.len() {
        let mut current = positions[start];
        let mut end = start + 1;
        while end < positions.len() && positions[end] < current + thr {
            current = positions[end];
            end += 1;
        }
        result.push((positions[start], positions[end - 1] + 1));
        start = end;
    }
    result
}

fn peak_detection(counts: &[usize], threshold: usize) -> Vec<(usize, usize)> {
    let (mut res, is_in_peak, start) = counts.iter().enumerate().fold(
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
    if is_in_peak {
        res.push((start, counts.len()));
    }
    res
}
