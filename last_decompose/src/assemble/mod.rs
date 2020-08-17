mod chunked_read;
pub mod correct_reads;
mod ditch_graph;
use super::Entry;
pub use chunked_read::ChunkedRead;
use std::collections::{HashMap, HashSet};

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

fn major_component(reads: &[ChunkedRead]) -> HashSet<u8> {
    let mut positions: HashMap<_, HashMap<_, u32>> = HashMap::new();
    let mut totals: HashMap<_, u32> = HashMap::new();
    // Register
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *totals.entry(node.window_position).or_default() += 1;
        }
        if let Some(cl) = read.label {
            for node in read.nodes.iter() {
                *positions
                    .entry(node.window_position)
                    .or_default()
                    .entry(cl)
                    .or_default() += 1;
            }
        }
    }
    let mut bg: HashMap<_, u32> = HashMap::new();
    let mut positions: Vec<_> = positions.into_iter().collect();
    positions.sort_by_key(|x| x.0);
    for (pos, counts) in positions {
        let total = *totals.get(&pos).unwrap();
        if let Some((&argmax, max)) = counts.iter().max_by_key(|x| x.1) {
            if total * 3 / 5 < *max && total > 70 {
                debug!("{}\t{}\t{}\t{}", pos, total, max, argmax);
                *bg.entry(argmax).or_default() += 1;
            }
        }
    }
    bg.iter()
        .filter(|&(_, &count)| count >= 1)
        .map(|(&cl, _)| cl)
        .collect()
}

type AssembleResult = (
    Vec<(String, Option<u8>)>,
    gfa::GFA,
    Vec<bio_utils::fasta::Record>,
);
pub fn assemble_reads(reads: &[ChunkedRead], k: usize, thr: usize) -> AssembleResult {
    let mut reads = reads.to_vec();
    correct_reads::correct_reads(&mut reads);
    // Determine SVs which are merged into the backgrounds.
    let backgrounds = major_component(&reads);
    debug!("backgrounds:{:?}", backgrounds);
    let background_cluster = backgrounds.iter().min().cloned();
    // Change labels and forbiddens.
    // All the backgrouds SVs are merged into one cluster.
    reads
        .iter_mut()
        .filter(|read| read.label.is_some())
        .filter(|read| backgrounds.contains(&read.label.unwrap()))
        .for_each(|read| {
            read.label = background_cluster;
        });
    reads
        .iter_mut()
        .filter(|r| !r.forbidden.is_empty())
        .for_each(|read| {
            if read.forbidden.iter().any(|f| backgrounds.contains(f)) {
                read.forbidden.retain(|f| !backgrounds.contains(f));
                read.forbidden.push(background_cluster.unwrap());
            }
        });
    {
        let count: Vec<_> = reads.iter().map(|n| n.nodes.len()).collect();
        let hist = histgram_viz::Histgram::new(&count);
        debug!("Read length:\n{}", hist.format(40, 20));

        let count: Vec<_> = reads.iter().map(|n| n.nodes.len()).collect();
        let hist = histgram_viz::Histgram::new(&count);
        debug!("Corrected Read length:\n{}", hist.format(40, 20));
    }
    let mut dbg = de_bruijn_graph::DeBruijnGraph::from(&reads, k);
    {
        let mut counts: HashMap<Vec<u64>, usize> = HashMap::new();
        for node in dbg.nodes.iter() {
            let kmer: Vec<_> = node.kmer.iter().map(|n| n.0).collect();
            *counts.entry(kmer).or_default() += node.occ;
        }
        let counts: Vec<_> = counts.values().copied().collect();
        let hist = histgram_viz::Histgram::new(&counts);
        debug!("Node occ:\n{}", hist.format(40, 20));
    }
    dbg.clean_up_auto();
    let count: Vec<_> = dbg.nodes.iter().map(|n| n.occ as usize).collect();
    let hist = histgram_viz::Histgram::new(&count);
    debug!("Hist of node occ:\n{}", hist.format(40, 20));
    let labels: Vec<Option<u8>> = reads.iter().map(|r| r.label).collect();
    let forbs: Vec<&[u8]> = reads.iter().map(|r| r.forbidden.as_slice()).collect();
    let label_map = dbg.expand_color(&reads, thr, &labels, &forbs, background_cluster);
    debug!("Merged Cluster:{:?}", label_map);
    let background_cluster = background_cluster.map(|l| label_map.get(&l).copied().unwrap_or(l));
    // Determine forbbidens labels on each components.
    let forbiddens_on_each_cluster: HashMap<_, Vec<u8>> = {
        let temp: HashMap<_, HashMap<_, u32>> = reads.iter().fold(HashMap::new(), |mut y, read| {
            let cluster = match read.label {
                Some(res) => Some(label_map.get(&res).cloned().unwrap_or(res)),
                None => dbg
                    .assign_read(read)
                    .or_else(|| dbg.assign_read_by_unit(read))
                    .map(|x| x as u8)
                    .or(background_cluster),
            };
            if let Some(cluster) = cluster {
                for &f in read.forbidden.iter() {
                    *y.entry(cluster).or_default().entry(f).or_default() += 1;
                }
            }
            y
        });
        temp.into_iter()
            .map(|(cl, counts)| {
                let forbidden: Vec<_> = counts.iter().filter(|x| x.1 > &50).map(|x| *x.0).collect();
                (cl, forbidden)
            })
            .collect()
    };
    for (cl, forbs) in forbiddens_on_each_cluster.iter() {
        debug!("Forbiddens on {}:{:?}", cl, forbs);
    }
    let assignments: Vec<_> = reads
        .iter()
        .map(|r| {
            let f = &forbiddens_on_each_cluster;
            let id = r.id.clone();
            let cluster = match r.label {
                //Some(res) => Some(*label_map.get(&res).unwrap_or(&res)),
                Some(res) => match label_map.get(&res) {
                    // Map labels into cluster on de Bruijn graph.
                    Some(mapped) => {
                        // If the `res` is not allowed to integrate
                        // `mapped` cluster,
                        // We return the original label.
                        if f.contains_key(mapped) && f[mapped].contains(&res) {
                            Some(res)
                        } else {
                            Some(*mapped)
                        }
                    }
                    // If there's no corresponding component on the de Bruijn graph
                    None => Some(res),
                },
                // We predict the label of the reads...
                None => dbg
                    .assign_read(r)
                    .or_else(|| dbg.assign_read_by_unit(r))
                    .map(|x| x as u8)
                    .or(background_cluster),
            };
            (id, cluster)
        })
        .collect();
    // Rename assignments.
    let map: HashMap<_, _> = {
        let clusters: HashSet<_> = assignments.iter().filter_map(|&(_, x)| x).collect();
        let mut clusters: Vec<_> = clusters.into_iter().collect();
        clusters.sort();
        clusters
            .into_iter()
            .enumerate()
            .map(|(x, y)| (y, x as u8))
            .collect()
    };
    debug!("Map:{:?}", map);
    let assignments: Vec<_> = assignments
        .into_iter()
        .map(|(id, asn)| (id, asn.map(|a| map[&a])))
        .collect();
    {
        let labels: HashMap<_, HashSet<String>> = reads.iter().fold(HashMap::new(), |mut x, y| {
            if let Some(label) = y.label {
                x.entry(label).or_default().insert(y.id.clone());
            }
            x
        });
        let result: HashMap<_, HashSet<String>> =
            assignments.iter().fold(HashMap::new(), |mut x, (id, asn)| {
                if let Some(asn) = asn {
                    x.entry(asn).or_default().insert(id.clone());
                }
                x
            });
        for (asn, ids) in result.iter() {
            debug!("Cluster:{}\t{}", asn, ids.len());
        }
        for (init_cl, members) in labels {
            debug!("Init cluster:{}", init_cl);
            for (result_cl, ids) in result.iter() {
                let count = members.intersection(ids).count();
                if count > 2 {
                    debug!("\tRes:{}({} reads)", result_cl, count);
                }
            }
        }
    }
    let max_cluster = map.values().max().cloned().unwrap_or(0);
    assert_eq!(assignments.len(), reads.len());
    let header = gfa::Content::Header(gfa::Header::default());
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
    debug!("Assembling reads...");
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
    (assignments, gfa, contigs)
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
    debug!("{}-{}len is_circular={}", id, seq.len(), is_circular);
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
    graph.remove_tips();
    graph.collapse_buddle();
    // debug!("{}", graph);
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
