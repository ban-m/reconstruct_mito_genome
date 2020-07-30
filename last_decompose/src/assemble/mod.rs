mod chunked_read;
mod ditch_graph;
use super::Entry;
pub use chunked_read::ChunkedRead;

#[derive(Debug, Clone)]
pub struct DecomposedResult {
    pub assignments: Vec<(String, Option<u8>)>,
    pub gfa: gfa::GFA,
    pub contigs: Vec<bio_utils::fasta::Record>,
}

impl de_bruijn_graph::AsDeBruijnNode for chunked_read::Node {
    fn as_node(w: &[chunked_read::Node]) -> de_bruijn_graph::Node {
        let first = w.first().unwrap().get_tuple();
        let last = w.last().unwrap().get_tuple();
        let kmer: Vec<_> = if first < last {
            w.iter().map(chunked_read::Node::get_tuple).collect()
        } else {
            w.iter().rev().map(chunked_read::Node::get_tuple).collect()
        };
        de_bruijn_graph::Node::new(kmer)
    }
}
impl de_bruijn_graph::IntoDeBruijnNodes for chunked_read::ChunkedRead {
    fn into_de_bruijn_nodes(&self, k: usize) -> Vec<de_bruijn_graph::Node> {
        self.nodes
            .windows(k)
            .map(de_bruijn_graph::AsDeBruijnNode::as_node)
            .collect()
    }
}

pub fn assemble_reads(reads: &mut [ChunkedRead]) -> DecomposedResult {
    let mut dbg = de_bruijn_graph::DeBruijnGraph::from(reads, 3);
    dbg.resolve_bubbles(reads);
    let max_cluster = dbg.coloring();
    let header = gfa::Content::Header(gfa::Header::default());
    let assignments: Vec<_> = reads
        .iter()
        .map(|r| (r.id.clone(), dbg.assign_read(r).map(|x| x as u8)))
        .collect();
    assert_eq!(assignments.len(), reads.len());
    let mut header = vec![gfa::Record::from_contents(header, vec![])];
    let clusters: Vec<Vec<_>> = (0..max_cluster)
        .map(|cl| {
            reads
                .iter()
                .zip(assignments.iter())
                .filter_map(|(read, (id, asn))| {
                    assert_eq!(id, &read.id);
                    asn.map(|asn| (read, asn))
                })
                .filter(|&(_, asn)| asn == cl as u8)
                .map(|x| x.0)
                .collect()
        })
        .collect();
    let records: Vec<_> = clusters
        .iter()
        .enumerate()
        .map(|(cl, reads)| reads_to_gfa(reads, cl))
        .collect();
    let contigs: Vec<_> = clusters
        .iter()
        .enumerate()
        .filter_map(|(cl, reads)| reads_to_contig(cl, reads))
        .collect();
    header.extend(records.into_iter().flat_map(|e| e.1));
    let gfa = gfa::GFA::from_records(header);
    DecomposedResult {
        assignments,
        contigs,
        gfa,
    }
}

fn reads_to_contig(cl: usize, reads: &[&ChunkedRead]) -> Option<bio_utils::fasta::Record> {
    debug!("Constructing the {}-th ditch graph", cl);
    if reads.len() < 10 {
        debug!("Detected small group:{}", reads.len());
        debug!("The reads below are discarded.");
        debug!("ID:Length:Nodes");
        for read in reads.iter() {
            debug!("{}:{}", read.id, read.nodes.len());
        }
        return None;
    }
    let mut graph = ditch_graph::DitchGraph::new(&reads);
    graph.collapse_buddle();
    debug!("{}", graph);
    let (seq, is_circular) = graph.simple_path();
    let desc = if is_circular {
        Some("is_circular=true".to_string())
    } else {
        None
    };
    let id = format!("tig_{:04}", cl);
    Some(bio_utils::fasta::Record::with_data(&id, &desc, &seq))
}

fn reads_to_gfa(reads: &[&ChunkedRead], cl: usize) -> (usize, Vec<gfa::Record>) {
    debug!("Constructing the {}-th ditch graph", cl);
    if reads.len() < 10 {
        debug!("Detected small group:{}", reads.len());
        debug!("The reads below are discarded.");
        debug!("ID:Length:Nodes");
        for read in reads.iter() {
            debug!("{}:{}", read.id, read.nodes.len());
        }
        return (cl, vec![]);
    }
    let mut graph = ditch_graph::DitchGraph::new(&reads);
    graph.collapse_buddle();
    debug!("{}", graph);
    let mut records = vec![];
    let (nodes, edges, group) = graph.spell(cl);
    let nodes = nodes
        .into_iter()
        .map(gfa::Content::Seg)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(nodes);
    let edges = edges
        .into_iter()
        .map(gfa::Content::Edge)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(edges);
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    records.push(group);
    (cl, records)
}
