use std::collections::HashMap;
#[derive(Clone)]
pub struct DitchGraph<'a, 'b> {
    nodes: Vec<DitchNode<'a, 'b>>,
    index: HashMap<NodeIndex, usize>,
}

impl<'a, 'b> std::fmt::Display for DitchGraph<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        use histgram_viz::Histgram;
        let nodes = self.nodes.len();
        let edges = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        writeln!(f, "Node:{}, Edges:{}", nodes, edges)?;
        let occs: Vec<_> = self.nodes.iter().map(|n| n.nodes.len()).collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Node Occs:{}", hist.format(20, 20))?;
        let occs: Vec<_> = self
            .nodes
            .iter()
            .flat_map(|n| n.edges.iter().map(|e| e.edges.len()))
            .collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Edge Occs:{}", hist.format(20, 20))?;
        let degrees = {
            let mut degs: HashMap<usize, usize> = HashMap::new();
            for node in self.nodes.iter() {
                let degree = node.edges.len();
                *degs.entry(degree).or_default() += 1;
            }
            let mut degs = degs.into_iter().collect::<Vec<_>>();
            degs.sort_by_key(|x| x.0);
            degs.into_iter()
                .map(|(deg, count)| format!("{}:{}", deg, count))
                .collect::<Vec<_>>()
        };
        write!(f, "[{}]", degrees.join(","))
    }
}

impl<'a, 'b> std::fmt::Debug for DitchGraph<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edge = self.nodes.iter().map(|e| e.edges.len()).sum::<usize>();
        writeln!(f, "Nodes:{}\tEdges:{}\n", self.nodes.len(), edge)?;
        let lines: Vec<_> = self
            .nodes
            .iter()
            .map(|node| {
                let idx = self
                    .index
                    .get(&NodeIndex::new(node.window_position, node.cluster))
                    .unwrap();
                format!("{}-{}", idx, node)
            })
            .collect();
        write!(f, "{}", lines.join("\n"))
    }
}

/// Ditch node.
/// It never allocates, or explicitly copy the nodes inside it,
/// rather, it borrows the reference to the `definitions::Node`s.
/// Note that even if a read contains (x,y,h) -> (z,w,t) connection,
/// the node of (z, w, t) contains an edge toward (x,y,h)-node.
/// Here, x and z are unit, y and w are cluster, and h and t are
/// head or tail.
#[derive(Debug, Clone)]
pub struct DitchNode<'a, 'b> {
    window_position: u64,
    cluster: u64,
    // CAUTION!!!!!! The `unit` and `cluster` members should not be look-upped!
    nodes: Vec<&'a super::chunked_read::Node>,
    edges: Vec<DitchEdge<'b>>,
}

impl<'a, 'b> std::fmt::Display for DitchNode<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(
            f,
            "---------{}:{}---------",
            self.window_position, self.cluster
        )?;
        writeln!(f, "Seq:")?;
        for (idx, n) in self.nodes.iter().enumerate() {
            writeln!(f, "{:3}:{}", idx, n.seq)?;
        }
        let lines: Vec<_> = self.edges.iter().map(|e| format!("{}", e)).collect();
        write!(f, "Edges:\n{}", lines.join("\n"))
    }
}

impl<'a, 'b> DitchNode<'a, 'b> {
    fn new(window_position: u64, cluster: u64) -> Self {
        Self {
            window_position,
            cluster,
            nodes: vec![],
            edges: vec![],
        }
    }
}

#[derive(Debug, Clone)]
pub struct DitchEdge<'a> {
    from: usize,
    to: usize,
    from_position: Position,
    to_position: Position,
    // If true, the edge is `forward` one.
    // In other words, if true,
    // you can spell this edge by its label,
    // otherwise you should take rev-cmp of the label.
    edges: Vec<(&'a str, bool)>,
}

impl<'a> std::fmt::Display for DitchEdge<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}:{})->", self.from, self.from_position)?;
        write!(f, "({}:{})", self.to, self.to_position)?;
        let edgs: Vec<_> = self
            .edges
            .iter()
            .map(|x| format!("{}", x.0.len()))
            .collect();
        write!(f, " {} threads[{}].", edgs.len(), edgs.join(","))
    }
}

impl<'a> std::cmp::PartialEq for DitchEdge<'a> {
    fn eq<'b>(&self, other: &DitchEdge<'b>) -> bool {
        self.from == other.from
            && self.to == other.to
            && self.from_position == other.from_position
            && self.to_position == other.to_position
    }
}

impl<'a> std::cmp::Eq for DitchEdge<'a> {}

impl<'a> DitchEdge<'a> {
    fn new(from: Position, from_i: usize, to: Position, to_i: usize) -> Self {
        Self {
            from: from_i,
            to: to_i,
            from_position: from,
            to_position: to,
            edges: vec![],
        }
    }
    fn push(&mut self, x: (&'a str, bool)) {
        self.edges.push(x);
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub enum Position {
    Head,
    Tail,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = match self {
            Position::Head => 'H',
            Position::Tail => 'T',
        };
        write!(f, "{}", x)
    }
}

impl std::ops::Not for Position {
    type Output = Self;
    fn not(self) -> Self {
        match self {
            Position::Head => Position::Tail,
            Position::Tail => Position::Head,
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash)]
pub struct NodeIndex {
    window_position: u64,
    cluster: u64,
}
impl NodeIndex {
    pub fn new(window_position: u64, cluster: u64) -> Self {
        Self {
            window_position,
            cluster,
        }
    }
}

fn revcmp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
}

// I use this enumeration to
// tag a contig to nodes of a ditch graph.
// Sometimes, a node of a ditch graph
// has both end and start of the contig to
// each tip (head and tail). Thus,
// I explicitly express the position of each tags.
// A contig named `String`, the length of which is `usize`,
// ends at Position, or starts at `Position`, or both.
#[derive(Debug, Clone, Eq, PartialEq)]
enum ContigTag {
    Start(String, Position, usize),
    End(String, Position, usize),
    // Start position, end position.
    Both(String, Position, Position, usize),
    None,
}

impl<'a> DitchGraph<'a, 'a> {
    pub fn new(reads: &[&'a super::ChunkedRead]) -> Self {
        // Allocate all nodes.
        let (index, mut nodes) = {
            let mut index = HashMap::new();
            let mut nodes = vec![];
            for node in reads.iter().flat_map(|r| r.nodes.iter()) {
                let window_position = node.window_position as u64;
                let node_index = NodeIndex::new(window_position, node.cluster as u64);
                index.entry(node_index).or_insert_with(|| {
                    nodes.push(DitchNode::new(window_position, node.cluster as u64));
                    nodes.len() - 1
                });
                let node_index = *index.get(&node_index).unwrap();
                nodes[node_index].nodes.push(node);
            }
            (index, nodes)
        };
        // Append edges.
        for read in reads.iter() {
            if read.nodes.len() != read.edges.len() + 1 {
                debug!("{}\t{}\t{}", read.id, read.nodes.len(), read.edges.len());
            }
            for (pairs, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
                let (from, to) = (&pairs[0], &pairs[1]);
                let from_position = from.window_position as u64;
                let from_index = *index
                    .get(&NodeIndex::new(from_position, from.cluster as u64))
                    .unwrap();
                let to_position = to.window_position as u64;
                let to_index = *index
                    .get(&NodeIndex::new(to_position, to.cluster as u64))
                    .unwrap();
                let from_pos = if from.is_forward {
                    Position::Tail
                } else {
                    Position::Head
                };
                let to_pos = if to.is_forward {
                    Position::Head
                } else {
                    Position::Tail
                };
                let mut dedg = DitchEdge::new(from_pos, from_index, to_pos, to_index);
                match nodes[from_index].edges.iter_mut().find(|x| x == &&dedg) {
                    Some(res) => res.push((edge, true)),
                    None => {
                        dedg.push((edge, true));
                        nodes[from_index].edges.push(dedg);
                    }
                }
                let mut dedg = DitchEdge::new(to_pos, to_index, from_pos, from_index);
                match nodes[to_index].edges.iter_mut().find(|x| x == &&dedg) {
                    Some(res) => res.push((edge, false)),
                    None => {
                        dedg.push((edge, false));
                        nodes[to_index].edges.push(dedg);
                    }
                }
            }
        }
        Self { index, nodes }
    }
    fn enumerate_candidates(&self) -> Vec<(usize, Position)> {
        let mut selected = vec![false; self.nodes.len()];
        let mut primary_candidates: Vec<_> = (0..self.nodes.len())
            .filter_map(|i| {
                for &position in &[Position::Head, Position::Tail] {
                    if self.nodes[i]
                        .edges
                        .iter()
                        .filter(|e| e.from_position == position)
                        .count()
                        != 1
                    {
                        selected[i] = true;
                        return Some((i, position));
                    }
                }
                None
            })
            .collect();
        let secondary_candidates =
            (0..self.nodes.len())
                .filter(|&i| !selected[i])
                .filter_map(|i| {
                    for &position in &[Position::Head, Position::Tail] {
                        // The unwrap() never panics.
                        let grand_child = self.nodes[i]
                            .edges
                            .iter()
                            .find(|e| e.from_position == position)
                            .map(|e| {
                                let count = self.nodes[e.to]
                                    .edges
                                    .iter()
                                    .filter(|f| f.from_position == e.to_position)
                                    .count();
                                if count == 0 {
                                    debug!("{},{}", i, position);
                                    debug!("Edge:{}", e);
                                    for edge in &self.nodes[e.to].edges {
                                        debug!("Child:{}", edge)
                                    }
                                    panic!()
                                }
                                count
                            })
                            .unwrap();
                        if grand_child > 1 {
                            return Some((i, position));
                        }
                    }
                    None
                });
        primary_candidates.extend(secondary_candidates);
        primary_candidates
    }
    // Assemble the ditch graph.
    // In other words, it reduce the simple path, currently.
    pub fn spell(&self, cl: usize) -> (Vec<gfa::Segment>, Vec<gfa::Edge>, gfa::Group) {
        let mut arrived = vec![false; self.nodes.len()];
        let mut sids = vec![ContigTag::None; self.nodes.len()];
        let (mut g_segs, mut g_edges) = (vec![], vec![]);
        let candidates = self.enumerate_candidates();
        for (i, p) in candidates {
            if arrived[i] {
                continue;
            }
            debug!("First loop:{}", i);
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges) = self.traverse_from(&mut arrived, &mut sids, i, p, name);
            g_segs.push(contig);
            g_edges.extend(edges);
        }
        for i in 0..self.nodes.len() {
            if arrived[i] {
                continue;
            }
            debug!("Second loop:{}", i);
            let p = Position::Head;
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges) = self.traverse_from(&mut arrived, &mut sids, i, p, name);
            g_segs.push(contig);
            g_edges.extend(edges);
        }
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.eid.clone()))
            .collect();
        let uid = Some(format!("group-{}", cl));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group)
    }
    // Removing nodes and re-map all of the indices.
    fn remove_nodes(&mut self, to_remove: &[bool]) {
        let mapping = {
            let mut mapping = vec![];
            let mut index = 0;
            for &b in to_remove {
                mapping.push(index);
                index += !b as usize;
            }
            mapping
        };
        self.index.iter_mut().for_each(|(_, v)| *v = mapping[*v]);
        {
            let mut idx = 0;
            self.nodes.retain(|_| {
                idx += 1;
                !to_remove[idx - 1]
            });
        }
        self.nodes.iter_mut().for_each(|node| {
            node.edges.iter_mut().for_each(|e| {
                e.from = mapping[e.from];
                e.to = mapping[e.to];
            })
        });
    }
    pub fn collapse_buddle(&mut self) {
        let mut to_remove = vec![false; self.nodes.len()];
        for i in 0..self.nodes.len() {
            for &position in &[Position::Head, Position::Tail] {
                let edges: Vec<_> = self.nodes[i]
                    .edges
                    .iter()
                    .filter(|e| e.from_position == position && e.from == i)
                    .collect();
                if edges.len() > 1 {
                    // Never panic.
                    let pos = edges[0].to_position;
                    let unit = self.nodes[edges[0].to].window_position;
                    if edges
                        .iter()
                        .all(|n| n.to_position == pos && unit == self.nodes[n.to].window_position)
                    {
                        // The head side of the i-th node has two or mode cluster,
                        // which has the same distination.
                        // Check the collaption criteria, and collapse if possible.
                        eprintln!(
                            "Removing a bubble starting from {}-{}-{:?}",
                            self.nodes[i].window_position, self.nodes[i].cluster, position
                        );
                        for i in self.collapse_bubble_from(i, position) {
                            // eprintln!("Removing {}", i);
                            to_remove[i] = true;
                        }
                    }
                }
            }
        }
        self.remove_nodes(&to_remove);
    }
    fn collapse_bubble_from(&mut self, i: usize, position: Position) -> Vec<usize> {
        // Check the collapsing condition.
        let edges: Vec<_> = self.nodes[i]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .inspect(|e| assert!(e.from == i))
            .collect();
        // Check.
        assert!(edges.len() > 1);
        let (first_id, first_pos, first_unit) = {
            let first = edges.first().unwrap();
            let first_unit = self.nodes[first.to].window_position;
            (first.to, first.to_position, first_unit)
        };
        let merged_nodes_ids: Vec<_> = edges.iter().map(|e| e.to).collect();
        assert!(edges
            .iter()
            .all(|e| e.to_position == first_pos && first_unit == self.nodes[e.to].window_position));
        // Get the candidate child.
        let pos = !first_pos;
        let (child_id, child_pos) = match self.nodes[first_id]
            .edges
            .iter()
            .find(|e| e.from_position == pos && e.from == first_id)
        {
            Some(child) => (child.to, child.to_position),
            None => return vec![],
        };
        // Check all distinations of `i` has only one parent, `i`.
        // Check all otherside of the distinations of `i` has only one child.
        let is_collapsable_bubble = edges.iter().all(|e| {
            let res = self.nodes[e.to].edges.iter().all(|f| {
                (f.to_position == position && f.to == i)
                    || (f.to_position == child_pos && f.to == child_id)
            });
            res && self.nodes[e.to].edges.len() == 2
        });
        if !is_collapsable_bubble {
            vec![]
        } else {
            // Add to the `first_id` from other edges.
            let remove_nodes: Vec<_> = edges.into_iter().skip(1).map(|e| e.to).collect();
            let mut node_result = vec![];
            let mut edge_result: Vec<(&str, bool)> = vec![];
            for &remove_node in remove_nodes.iter() {
                node_result.append(&mut self.nodes[remove_node].nodes);
                self.nodes[remove_node]
                    .edges
                    .iter_mut()
                    .for_each(|ditch_edge| {
                        edge_result.append(&mut ditch_edge.edges);
                    });
            }
            self.nodes[first_id].nodes.extend(node_result);
            self.nodes[first_id]
                .edges
                .iter_mut()
                .find(|ditch_edge| ditch_edge.to == i && ditch_edge.to_position == position)
                .unwrap()
                .edges
                .extend(edge_result);
            // Change edges from nodes[i]
            let edge_result = self.nodes[i]
                .edges
                .iter_mut()
                .filter(|e| e.from_position == position)
                .skip(1)
                .fold(vec![], |mut x, ditch_edge| {
                    x.append(&mut ditch_edge.edges);
                    x
                });
            self.nodes[i]
                .edges
                .iter_mut()
                .find(|e| e.from_position == position)
                .unwrap()
                .edges
                .extend(edge_result);
            self.nodes[i]
                .edges
                .retain(|ditch_edge| !ditch_edge.edges.is_empty());
            let edge_result = self.nodes[child_id]
                .edges
                .iter_mut()
                .filter(|e| merged_nodes_ids.contains(&e.to))
                .skip(1)
                .fold(vec![], |mut x, d_edge| {
                    x.append(&mut d_edge.edges);
                    x
                });
            // Never panic.
            self.nodes[child_id]
                .edges
                .iter_mut()
                .find(|e| merged_nodes_ids.contains(&e.to))
                .unwrap()
                .edges
                .extend(edge_result);
            self.nodes[child_id]
                .edges
                .retain(|ditch_edge| !ditch_edge.edges.is_empty());
            remove_nodes
        }
    }
    // Traverse from the given `start` node.
    fn traverse_from(
        &self,
        arrived: &mut [bool],
        sids: &mut [ContigTag],
        start: usize,
        start_position: Position,
        seqname: String,
    ) -> (gfa::Segment, Vec<gfa::Edge>) {
        // Find edges.
        let mut edges: Vec<_> = self.nodes[start]
            .edges
            .iter()
            .filter(|e| e.from_position == start_position)
            .filter_map(|e| {
                assert!(e.from == start);
                let (sid2, beg2) = match sids[e.to] {
                    ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(0, false);
                let a = None;
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            })
            .collect();
        let (mut node, mut position) = (start, start_position);
        let mut seq = String::new();
        // Start traveresing.
        let mut unit_names = vec![];
        loop {
            unit_names.push((self.nodes[node].window_position, self.nodes[node].cluster));
            arrived[node] = true;
            // Move forward.
            match position {
                Position::Head => seq += &self.nodes[node].nodes[0].seq,
                Position::Tail => {
                    let s = &self.nodes[node].nodes[0].seq;
                    seq += &revcmp(s);
                }
            };
            position = !position;
            // Check.
            let num_edges = self.nodes[node]
                .edges
                .iter()
                .filter(|e| e.from_position == position)
                .count();
            if num_edges == 0 || num_edges > 1 {
                break;
            }
            // There is only one child.
            let selected_edge = self.nodes[node]
                .edges
                .iter()
                .find(|e| e.from_position == position)
                .unwrap();
            assert_eq!(selected_edge.from, node);
            assert_eq!(selected_edge.from_position, position);
            // Succeed along with the edge.
            // if let Some(&(ref edge, is_forward)) =
            let &(ref edge, is_forward) = selected_edge.edges.first().unwrap();
            if is_forward {
                seq += &edge;
            } else {
                seq += &revcmp(&edge);
            }
            let (next, next_position) = (selected_edge.to, selected_edge.to_position);
            // Check the number of child.
            let num_children = self.nodes[next]
                .edges
                .iter()
                .filter(|e| e.from_position == next_position)
                .count();
            if num_children >= 2 || arrived[next] {
                break;
            }
            // Jump to the next node.
            assert!(num_children == 1);
            node = next;
            position = next_position;
        }
        let unit_names: Vec<_> = unit_names
            .iter()
            .map(|(u, c)| format!("{}:{}", u, c))
            .collect();
        debug!("{}\t{}", seqname, unit_names.join("\t"));
        // Register start and tail node.
        if start == node {
            sids[node] = ContigTag::Both(seqname.clone(), start_position, position, seq.len());
        } else {
            sids[start] = ContigTag::Start(seqname.clone(), start_position, seq.len());
            sids[node] = ContigTag::End(seqname.clone(), position, seq.len());
        }
        let seg = gfa::Segment::from(seqname.clone(), seq.len(), Some(seq.clone()));
        // Add gfa edges.
        let tail_edges = self.nodes[node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .filter_map(|e| {
                assert!(e.from == node);
                let (sid2, beg2) = match sids[e.to] {
                    ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(seq.len(), true);
                let a = None;
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            });
        edges.extend(tail_edges);
        (seg, edges)
    }
    pub fn simple_path(&self) -> (Vec<u8>, bool) {
        let (start, position, max, is_circular) = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, n)| {
                let tail_num = n
                    .edges
                    .iter()
                    .filter(|e| e.from_position == Position::Tail)
                    .count();
                let head_num = n
                    .edges
                    .iter()
                    .filter(|e| e.from_position == Position::Head)
                    .count();
                if tail_num == 0 {
                    Some((idx, Position::Tail, n.nodes.len(), false))
                } else if head_num == 0 {
                    Some((idx, Position::Head, n.nodes.len(), false))
                } else {
                    None
                }
            })
            .max_by_key(|n| n.2)
            .unwrap_or_else(|| {
                let (start, max) = self
                    .nodes
                    .iter()
                    .enumerate()
                    .max_by_key(|n| n.1.nodes.len())
                    .unwrap();
                (start, Position::Head, max.nodes.len(), true)
            });
        debug!("{}-{}(Occ:{})", start, position, max);
        let mut node = start;
        let mut position = position;
        let mut seq = String::new();
        let mut arrived = vec![false; self.nodes.len()];
        loop {
            arrived[node] = true;
            match position {
                Position::Head => seq += &self.nodes[node].nodes[0].seq,
                Position::Tail => {
                    let s = &self.nodes[node].nodes[0].seq;
                    seq += &revcmp(s);
                }
            };
            position = !position;
            let selected_edge = self.nodes[node]
                .edges
                .iter()
                .filter(|e| e.from_position == position)
                .filter(|e| !arrived[e.to])
                .max_by_key(|e| self.nodes[e.to].nodes.len());
            let selected_edge = match selected_edge {
                Some(edge) => edge,
                None => break,
            };
            assert_eq!(selected_edge.from, node);
            assert_eq!(selected_edge.from_position, position);
            let &(ref edge, is_forward) = selected_edge.edges.first().unwrap();
            if is_forward {
                seq += &edge;
            } else {
                seq += &revcmp(&edge);
            }
            node = selected_edge.to;
            position = selected_edge.to_position;
        }
        (seq.as_bytes().to_vec(), is_circular)
    }
}
