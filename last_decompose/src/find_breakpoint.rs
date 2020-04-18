use super::ERead;
use last_tiling::repeat::RepeatPairs;
use last_tiling::{Contigs, LastTAB};
use serde::{Deserialize, Serialize};
use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::HashSet;
use std::fmt::{Debug, Display, Formatter};
//const MERGE_THR: usize = 500;
const MERGE_THR: usize = 5;
pub const COVERAGE_THR: usize = 25;
const WINDOW: usize = 3;
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cluster {
    pub id: usize,
    pub members: Vec<Member>,
}

impl Cluster {
    pub fn ids(&self) -> HashSet<String> {
        self.members.iter().fold(HashSet::new(), |mut x, y| {
            x.extend(y.cr.reads().iter().cloned());
            x
        })
    }
    pub fn ranges(&self) -> Vec<(u16, (i32, i32))> {
        self.members
            .iter()
            .flat_map(|mem| match mem.cr {
                CriticalRegion::CP(ref cp) => vec![
                    (cp.contig1().contig(), cp.contig1().range()),
                    (cp.contig2().contig(), cp.contig2().range()),
                ],
                CriticalRegion::CR(ref cr) => vec![(cr.contig().contig(), cr.contig().range())],
            })
            .collect()
    }
    pub fn overlap(&self, range: (u16, u16, u16)) -> bool {
        self.members.iter().any(|m| m.cr.overlap(range))
    }
}

impl ReadClassify for Cluster {
    fn is_spanned_by(&self, r: &ERead) -> bool {
        self.members.iter().any(|m| m.cr.is_spanned_by(r))
    }
    fn contains(&self, r: &ERead) -> bool {
        self.members.iter().any(|m| m.cr.contains(r))
    }
    fn along_with(&self, r: &ERead) -> bool {
        self.members.iter().any(|m| m.cr.along_with(r))
    }
    fn has(&self, id: &str) -> bool {
        self.members.iter().any(|m| m.cr.has(id))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Member {
    pub cr: CriticalRegion,
    pub cluster: usize,
}

fn get_forbids_cluster<'a>(
    reads: &'a [ERead],
    clusters: &[Vec<CriticalRegion>],
) -> HashMap<&'a str, Vec<usize>> {
    reads
        .iter()
        .map(|r| {
            let forbs: Vec<_> = clusters
                .iter()
                .enumerate()
                .filter(|(_, clsuters)| clsuters.iter().any(|c| c.is_spanned_by(r)))
                .map(|(id, _)| id)
                .collect();
            (r.id(), forbs)
        })
        .collect()
}

fn intersection_of(
    crs: &[Vec<CriticalRegion>],
    (i, j): (usize, usize),
    forbids: &HashMap<&str, Vec<usize>>,
) -> usize {
    let xs: HashSet<_> = crs[i]
        .iter()
        .flat_map(|cr| cr.reads())
        .filter(|&id| !forbids[id.as_str()].contains(&j))
        .collect();
    let ys: HashSet<_> = crs[j]
        .iter()
        .flat_map(|cr| cr.reads())
        .filter(|&id| !forbids[id.as_str()].contains(&i))
        .collect();
    xs.intersection(&ys).count()
}

fn merge(
    crs: &mut Vec<Vec<CriticalRegion>>,
    (i, j): (usize, usize),
    forbiddens: &HashMap<&str, Vec<usize>>,
) {
    let mut j_clusters = crs.remove(j);
    for j_cluster in j_clusters.iter_mut() {
        let reads: HashSet<String> = j_cluster
            .reads()
            .iter()
            .filter(|id| !forbiddens[id.as_str()].contains(&i))
            .cloned()
            .collect();
        j_cluster.replace_reads(reads);
    }
    for i_cluster in crs[i].iter_mut() {
        let reads: HashSet<String> = i_cluster
            .reads()
            .iter()
            .filter(|id| !forbiddens[id.as_str()].contains(&j))
            .cloned()
            .collect();
        i_cluster.replace_reads(reads);
    }
    crs[i].extend(j_clusters);
}

pub fn initial_clusters(
    reads: &[ERead],
    contigs: &Contigs,
    repeats: &[RepeatPairs],
    alignments: &[LastTAB],
) -> Vec<Cluster> {
    let crs: Vec<_> = critical_regions(reads, contigs, repeats, alignments);
    let mut crs: Vec<_> = crs.into_iter().map(|e| vec![e]).collect();
    'merge: loop {
        let len = crs.len();
        let forbiddens = get_forbids_cluster(reads, &crs);
        debug!("Current Cluster:{}", len);
        for i in 0..len {
            for j in (i + 1)..len {
                let intersection = intersection_of(&crs, (i, j), &forbiddens);
                if intersection > COVERAGE_THR {
                    debug!("Cluster {} and {} share {} reads.", i, j, intersection);
                    merge(&mut crs, (i, j), &forbiddens);
                    continue 'merge;
                }
            }
        }
        break;
    }
    debug!("Resulting in {} clusters.", crs.len());
    crs.into_iter()
        .enumerate()
        .map(|(id, members)| {
            let members = members
                .into_iter()
                .map(|cr| Member { cr, cluster: id })
                .collect();
            Cluster { id, members }
        })
        .collect()
}

/// How many reads are needed to construct a separate class.
/// Note that it should be more than 1, since there is a
/// danger of chimeric reads.
const CHECK_THR: u16 = 4;
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
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
    fn has(&self, id: &str) -> bool;
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
    fn has(&self, r: &str) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.has(r),
            CriticalRegion::CR(ref cr) => cr.has(r),
        }
    }
}
impl CriticalRegion {
    pub fn reads(&self) -> &HashSet<String> {
        match self {
            CriticalRegion::CP(ref cp) => &cp.reads,
            CriticalRegion::CR(ref cr) => &cr.reads,
        }
    }
    pub fn replace_reads(&mut self, reads: HashSet<String>) {
        match self {
            CriticalRegion::CP(cp) => cp.reads = reads,
            CriticalRegion::CR(cr) => cr.reads = reads,
        }
    }
    fn overlap(&self, range: (u16, u16, u16)) -> bool {
        match self {
            CriticalRegion::CP(ref cp) => cp.overlap(range),
            CriticalRegion::CR(ref cr) => cr.overlap(range),
        }
    }
    pub fn contig_pair(&self) -> Option<&ContigPair> {
        match self {
            CriticalRegion::CP(ref cp) => Some(cp),
            CriticalRegion::CR(_) => None,
        }
    }
    pub fn confluent_region(&self) -> Option<&ConfluentRegion> {
        match self {
            CriticalRegion::CP(_) => None,
            CriticalRegion::CR(ref cr) => Some(cr),
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
    fn overlap(&self, (contig, start, end): (u16, u16, u16)) -> bool {
        let overlap = !(end as u16 <= self.start_unit || self.end_unit < start as u16);
        contig == self.contig && overlap
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
        let s = self.start_unit.max(CHECK_THR) - CHECK_THR;
        let t = self.end_unit + CHECK_THR;
        let (s_thr, t_thr) = (s.max(CHECK_THR) - CHECK_THR, t + CHECK_THR);
        if t_thr < self.longest {
            r.does_touch(self.contig, s_thr, s) && r.does_touch(self.contig, t, t_thr)
        } else {
            r.does_touch(self.contig, s_thr, s) && r.does_touch(self.contig, 0, CHECK_THR)
        }
    }
}

#[derive(Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct ContigPair {
    contig1: Position,
    contig2: Position,
    reads: HashSet<String>,
}

impl Debug for ContigPair {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        writeln!(f, "contig1:{:?}", self.contig1)?;
        writeln!(f, "contig2:{:?}", self.contig2)?;
        writeln!(f, "Number of reads:{}", self.reads.len())
    }
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
        self.reads.contains(r.id())
    }
    fn has(&self, id: &str) -> bool {
        self.reads.contains(id)
    }
}

impl ContigPair {
    pub fn new(contig1: Position, contig2: Position, reads: HashSet<String>) -> Self {
        Self {
            contig1,
            contig2,
            reads,
        }
    }
    pub fn contig1(&self) -> &Position {
        &self.contig1
    }
    pub fn contig2(&self) -> &Position {
        &self.contig2
    }
    fn overlap(&self, range: (u16, u16, u16)) -> bool {
        self.contig1.overlap(range) || self.contig2.overlap(range)
    }
    fn overlap_with(&self, cr: &ConfluentRegion) -> bool {
        let Position {
            contig,
            start_unit,
            end_unit,
            ..
        } = cr.pos;
        self.overlap((contig, start_unit - CHECK_THR, end_unit + CHECK_THR))
    }
    pub fn reads(&self) -> &HashSet<String> {
        &self.reads
    }
}

#[derive(Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct ConfluentRegion {
    // Contains start and end position.
    pos: Position,
    reads: HashSet<String>,
}

impl Debug for ConfluentRegion {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        writeln!(f, "pos:{:?}", self.pos)?;
        writeln!(f, "Number of reads:{}", self.reads.len())
    }
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
    fn new(pos: Position, reads: HashSet<String>) -> Self {
        Self { pos, reads }
    }
    pub fn contig(&self) -> &Position {
        &self.pos
    }
    fn overlap(&self, range: (u16, u16, u16)) -> bool {
        self.pos.overlap(range)
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
        self.reads.contains(r.id())
    }
    fn has(&self, id: &str) -> bool {
        self.reads.contains(id)
    }
}

/// Return critical regions.
pub fn critical_regions(
    reads: &[ERead],
    contigs: &Contigs,
    _repeats: &[RepeatPairs],
    _alignments: &[LastTAB],
) -> Vec<CriticalRegion> {
    let contig_pairs = contigpair_position(reads, contigs);
    let confluent_regions = confluent_position(reads, contigs, last_tiling::UNIT_SIZE);
    let confluent_regions: Vec<_> = confluent_regions
        .into_iter()
        .filter(|cr| contig_pairs.iter().all(|cp| !cp.overlap_with(cr)))
        .map(CriticalRegion::CR)
        .collect();
    contig_pairs
        .into_iter()
        .map(CriticalRegion::CP)
        .chain(confluent_regions)
        .collect()
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
            let former_from_up = w[0].unit < w[1].unit;
            // If w[2] > w[3], it is to upstream.
            let latter_from_up = w[2].unit > w[3].unit;
            match (start < end, former_from_up, latter_from_up) {
                (true, true, x) => from_up_counts[c1][start].push((c2 as u16, end, x)),
                (true, false, x) => to_down_counts[c1][start].push((c2 as u16, end, x)),
                (false, x, true) => from_up_counts[c2][end].push((c1 as u16, start, x)),
                (false, x, false) => to_down_counts[c2][end].push((c1 as u16, start, x)),
            }
        }
    }
    let mut contigpairs = vec![];
    let pairs = peak_call_contigpair(reads, from_up_counts, true, contigs);
    contigpairs.extend(pairs);
    let pairs = peak_call_contigpair(reads, to_down_counts, false, contigs);
    contigpairs.extend(pairs);
    contigpairs
}

// This function is broken!!!!
// FIXME!!!!!
fn peak_call_contigpair(
    reads: &[ERead],
    mut counts: Vec<Vec<Vec<(u16, usize, bool)>>>,
    from_upstream: bool,
    contigs: &Contigs,
) -> Vec<ContigPair> {
    let mut contigpairs = vec![];
    for (contig, jumps) in counts.iter_mut().enumerate() {
        let contig = contig as u16;
        let from_max = contigs.get_last_unit(contig).unwrap();
        let mut start_idx = 0;
        while start_idx < jumps.len() {
            // Search a heviest edge from this contig.
            let (from, to_c, to, direction) = match search_jump_start(jumps, start_idx) {
                Some(res) => res,
                None => break,
            };
            // Determine the range of the start position ad the range of
            // the end position.
            let (from_start, from_end, to_start, to_end) =
                search_jump_area(jumps, from, to_c, to, direction);
            // Collect the edges from (start,end) to (start_dst, end_dst)
            let counts = collect_edges_in(
                jumps, from_start, from_end, to_c, to_start, to_end, direction,
            );
            // and remove them from counts.
            if counts > COVERAGE_THR {
                remove_edges(
                    jumps, from_start, from_end, to_c, to_start, to_end, direction,
                );
                let (belong_reads, total) = get_belong_reads(
                    reads,
                    (contig, from_start, from_end, from_upstream),
                    (to_c, to_start, to_end, direction),
                );
                assert_eq!(total, counts);
                use Direction::*;
                let from_direction = if from_upstream { UpStream } else { DownStream };
                let to_direction = if direction { UpStream } else { DownStream };
                let to_max = contigs.get_last_unit(to_c as u16).unwrap();
                let (contig, from_start, from_end) =
                    (contig as u16, from_start as u16, from_end as u16);
                let pos1 = Position::new(contig, from_start, from_end, from_direction, from_max);
                let (to_start, to_end) = (to_start as u16, to_end as u16);
                let pos2 = Position::new(to_c, to_start, to_end, to_direction, to_max);
                contigpairs.push(ContigPair::new(pos1, pos2, belong_reads));
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
        if let Some((&(x, y, z), _)) = counts
            .iter()
            .filter(|&(_, &count)| count > JUMP_THR)
            .min_by_key(|&(&(_, pos, _), _)| pos)
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
        (min, max + 1)
    };
    let (mut from_start, mut from_end) = (from, from + 1);
    let (mut to_start, mut to_end) = get_range(from_start, to, to + 1);
    const STEP_SIZE: usize = 6;
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

fn get_belong_reads(
    reads: &[ERead],
    (f_c, f_start, f_end, f_dir): (u16, usize, usize, bool),
    (t_c, t_start, t_end, t_dir): (u16, usize, usize, bool),
) -> (HashSet<String>, usize) {
    // debug!(
    //     "{}:{}-{}=>{}:{}-{}",
    //     f_c, f_start, f_end, t_c, t_start, t_end
    // );
    let mut total = 0;
    let is_along_with = |read: &&ERead| -> bool {
        let jumps = read.seq().windows(4).filter(is_jumping);
        let count = jumps
            .filter(|w| {
                let (start, end) = (w[1].unit as usize, w[2].unit as usize);
                let (c1, c2) = (w[1].contig, w[2].contig);
                let former_direction = w[0].unit < w[1].unit;
                let latter_direction = w[2].unit > w[3].unit;
                if start < end {
                    // former <=> from, latter <=> to
                    let direction = former_direction == f_dir && latter_direction == t_dir;
                    let contig = f_c == c1 && t_c == c2;
                    let range = f_start <= start && start < f_end && t_start <= end && end < t_end;
                    direction & contig & range
                } else {
                    // former <=> to, latter <=> from.
                    let direction = former_direction == t_dir && latter_direction == f_dir;
                    let contig = t_c == c1 && f_c == c2;
                    let range = f_start <= end && end < f_end && t_start <= start && start < t_end;
                    direction & contig & range
                }
            })
            .count();
        total += count;
        count != 0
    };
    let result: HashSet<_> = reads
        .iter()
        .filter(is_along_with)
        .map(|r| r.id().to_string())
        .collect();
    (result, total)
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

/// Enumerate confluent region in the references.
pub fn confluent_position(
    reads: &[ERead],
    contigs: &Contigs,
    _unit_size: usize,
) -> Vec<ConfluentRegion> {
    let mut start_stop_count: HashMap<u16, Vec<_>> = contigs
        .get_last_units()
        .into_iter()
        .enumerate()
        .map(|(id, len)| (id as u16, vec![(0, 0); len as usize + 1]))
        .collect();
    for read in reads.iter().filter(|r| r.seq().len() > 1) {
        if let Some((first, last)) = read.get_edges() {
            {
                let (first, is_downstream) = first;
                let res = start_stop_count.get_mut(&first.contig).unwrap();
                if is_downstream {
                    res[first.unit as usize].0 += 1;
                } else {
                    res[first.unit as usize].1 += 1;
                }
            }
            {
                let (last, is_downstream) = last;
                let res = start_stop_count.get_mut(&last.contig).unwrap();
                if is_downstream {
                    res[last.unit as usize].0 += 1;
                } else {
                    res[last.unit as usize].1 += 1;
                }
            }
        }
    }
    let mut result = vec![];
    for (id, counts) in start_stop_count {
        let downstream: Vec<_> = counts.iter().map(|x| x.0).collect();
        let chunks = peak_call(&downstream).into_iter().map(|(s, e)| {
            let count = downstream[s..e].iter().sum::<usize>();
            (id, s, e, true, count)
        });
        result.extend(chunks);
        let upstream: Vec<_> = counts.iter().map(|x| x.1).collect();
        let chunks = peak_call(&upstream).into_iter().map(|(s, e)| {
            let count = upstream[s..e].iter().sum::<usize>();
            (id, s, e, false, count)
        });
        result.extend(chunks);
    }
    result
        .iter()
        .filter_map(|&(id, start, end, to_downstream, _)| {
            let reads = get_confluent_reads(&reads, (id, start, end, to_downstream));
            let (start, end) = (start as u16, end as u16);
            let max = contigs.get_last_unit(id).unwrap();
            use Direction::*;
            let di = if to_downstream { DownStream } else { UpStream };
            let is_edge = (end < OFFSET as u16 && to_downstream)
                || (max < OFFSET as u16 + start && !to_downstream);
            if is_edge {
                None
            } else {
                let cr = ConfluentRegion::new(Position::new(id, start, end, di, max), reads);
                Some(cr)
            }
        })
        .collect()
}

fn get_confluent_reads(
    reads: &[ERead],
    (id, start, end, to_down): (u16, usize, usize, bool),
) -> HashSet<String> {
    let is_unit_match = |(unit, is_downstream): (&crate::CUnit, bool)| {
        unit.contig == id
            && start as u16 <= unit.unit
            && unit.unit < end as u16
            && to_down == is_downstream
    };
    let is_match = |r: &&ERead| match r.get_edges() {
        Some((first, last)) => is_unit_match(first) || is_unit_match(last),
        None => false,
    };
    reads
        .iter()
        .filter(is_match)
        .map(|x| x.id().to_string())
        .collect()
}

fn peak_call(counts: &[usize]) -> Vec<(usize, usize)> {
    let average = counts.iter().sum::<usize>() / counts.len() * 2;
    debug!("AVE:{}", average);
    let thr = average.max(COVERAGE_THR);
    let positions: Vec<_> = counts
        .iter()
        .enumerate()
        .filter(|(_, &count)| count >= thr)
        .map(|(idx, _)| idx)
        .collect();
    merge_neighbors(positions, MERGE_THR)
        .into_iter()
        .filter(|&(s, t)| counts[s..t].iter().sum::<usize>() / (t - s + 1) > thr)
        .collect()
}

fn merge_neighbors(positions: Vec<usize>, thr: usize) -> Vec<(usize, usize)> {
    let mut result: Vec<(usize, usize)> = vec![];
    let mut start = 0;
    while start < positions.len() {
        let mut current = positions[start];
        let mut end = start;
        while end < positions.len() && positions[end] < current + thr {
            current = positions[end];
            end += 1;
        }
        let start_pos = positions[start];
        let end_pos = positions[end - 1] + WINDOW;
        result.push((start_pos, end_pos));
        start = end;
    }
    result
}
