use last_decompose::d3_data::Summary;
use std::collections::HashSet;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::fs::File;
    use std::io::BufReader;
    let Summary {
        clusters, reads, ..
    } = serde_json::de::from_reader(BufReader::new(File::open(&args[2])?)).unwrap();
    let target: usize = args[1].parse().unwrap();
    let cluster = clusters.into_iter().find(|x| x.id == target).unwrap();
    println!("Target Cluster: {}", target);
    use last_decompose::d3_data::Unit;
    let meets_end = |confluent: &last_decompose::find_breakpoint::ConfluentRegion| {
        reads
            .iter()
            .filter(|read| read.cluster() as usize == target)
            .filter(|read| {
                let units = read.units().iter().filter_map(|unit| match unit {
                    Unit::G(_) => None,
                    Unit::E(contig, pos) => Some((contig, pos)),
                });
                let (start, end) = confluent.contig().range();
                let start = start as u16;
                let end = end as u16;
                let c_contig = confluent.contig().contig();
                let first_boundary = {
                    let (&contig, &pos) = units.clone().nth(0).unwrap();
                    c_contig == contig && start <= pos && pos < end
                };
                let last_boundary = {
                    let (&contig, &pos) = units.clone().last().unwrap();
                    c_contig == contig && start <= pos && pos < end
                };
                first_boundary || last_boundary
            })
            .map(|r| r.name().to_string())
            .collect::<HashSet<_>>()
    };
    let meets_boundaries: Vec<(HashSet<_>, _)> = cluster
        .members
        .iter()
        .filter_map(|cl| cl.cr.confluent_region())
        .map(|confluent| (meets_end(confluent), confluent))
        .collect();
    for (rs, confluent) in meets_boundaries.iter() {
        println!("{:?} has {} reads", confluent, rs.len());
    }
    let union: HashSet<_> = meets_boundaries
        .iter()
        .fold(HashSet::new(), |mut acc, (x, _)| {
            acc.extend(x.iter().cloned());
            acc
        });
    println!("Union:{}", union.len());
    let intersection = {
        let (acc, _) = meets_boundaries[0].clone();
        meets_boundaries[0..]
            .iter()
            .fold(acc, |acc, (x, _)| acc.intersection(x).cloned().collect())
    };
    println!("Intersection:{}", intersection.len());
    // for read in reads {
    //     if intersection.contains(read.name()) {
    //         let units: Vec<_> = read
    //             .units()
    //             .iter()
    //             .map(|u| match u {
    //                 Unit::E(x, y) => format!("{}:{}-", x, y),
    //                 Unit::G(x) => format!("{}-", x),
    //             })
    //             .collect();
    //         println!("{}:{}\n", read.name(), units.join(""));
    //     }
    // }
    Ok(())
}
