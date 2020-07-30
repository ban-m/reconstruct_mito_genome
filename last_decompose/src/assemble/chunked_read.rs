#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChunkedRead {
    pub id: String,
    pub desc: Option<String>,
    pub nodes: Vec<Node>,
    pub edges: Vec<String>,
    // Label detemined by intial SV detection.
    pub label: Option<u8>,
    // Forbidden clusters.
    pub forbidden: Vec<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node {
    pub seq: String,
    pub cluster: u8,
    pub window_position: usize,
    pub is_forward: bool,
}

impl ChunkedRead {
    pub fn from<'a>(
        r: &last_tiling::EncodedRead,
        label: Option<&u8>,
        forbidden: Option<&Vec<u8>>,
        entries: &[super::Entry<'a>],
    ) -> Self {
        // assert!(entries.is_sorted_by_key(|e| e.read_position));
        use last_tiling::unit::ChunkedUnit::{En, Gap};
        let node_positions: Vec<_> = entries
            .iter()
            .map(|entry| {
                let (ctg, start, end) = entry.window_range;
                let mut units = r
                    .seq()
                    .iter()
                    .enumerate()
                    .skip_while(|(_, unit)| match unit {
                        En(e) => !(ctg == e.contig && start <= e.unit && e.unit < end),
                        Gap(_) => true,
                    })
                    .take_while(|(_, unit)| match unit {
                        En(e) => (ctg == e.contig && start <= e.unit && e.unit < end),
                        Gap(_) => true,
                    });
                // Never panic!
                let start_idx = units.next().unwrap().0;
                let end_idx = units.last().map(|x| x.0).unwrap_or(start_idx + 1);
                (start_idx, end_idx)
            })
            .collect();
        let nodes: Vec<_> = node_positions
            .iter()
            .zip(entries.iter())
            .map(|(&(s, e), entry)| {
                let is_forward = r.is_forward_wrt(entry.window_range.0).unwrap();
                let seq: Vec<_> = r.seq()[s..e]
                    .iter()
                    .flat_map(|unit| match unit {
                        En(e) if e.is_forward => e.bases.as_bytes().to_vec(),
                        En(e) => bio_utils::revcmp(e.bases.as_bytes()),
                        Gap(g) => g.bases().to_vec(),
                    })
                    .collect();
                let seq = String::from_utf8_lossy(&seq);
                Node {
                    seq: seq.to_string(),
                    cluster: entry.assignment,
                    window_position: entry.window,
                    is_forward,
                }
            })
            .collect();
        let edges: Vec<_> = node_positions
            .windows(2)
            .map(|w| {
                let start = w[0].1;
                let end = w[1].0;
                if start < end {
                    let seq: Vec<_> = r.seq()[start..end]
                        .iter()
                        .flat_map(|unit| match unit {
                            En(e) if e.is_forward => e.bases.as_bytes().to_vec(),
                            En(e) => bio_utils::revcmp(e.bases.as_bytes()),
                            Gap(g) => g.bases().to_vec(),
                        })
                        .collect();
                    String::from_utf8_lossy(&seq).to_string()
                } else {
                    String::new()
                }
            })
            .collect();
        let id = r.id.clone();
        let desc = r.desc.clone();
        let label = label.copied();
        let forbidden = forbidden.cloned().unwrap_or_else(Vec::new);
        Self {
            id,
            desc,
            nodes,
            edges,
            forbidden,
            label,
        }
    }
}
impl Node {
    pub fn get_tuple(&self) -> (u64, u64) {
        (self.window_position as u64, self.cluster as u64)
    }
}
