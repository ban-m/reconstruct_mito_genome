use last_decompose::d3_data::convert_to_d3_data_with_assign;
use last_decompose::find_breakpoint::Cluster;
use last_tiling::{Contigs, EncodedRead};
use log::debug;
use std::collections::HashMap;
pub mod template;

pub fn last_alignment_train<P: AsRef<std::path::Path>>(
    query: &P,
    refr: &P,
    thread: usize,
) -> Option<Vec<last_tiling::LastTAB>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 10_000;
    let mut c_dir = std::env::current_dir().ok()?;
    use std::io::{BufWriter, Write};
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir).ok()?;
    let reference = refr.as_ref().to_str()?;
    let reads = query.as_ref().to_str()?;
    let db_name = {
        let mut temp = c_dir.clone();
        temp.push("reference");
        temp.into_os_string().into_string().unwrap()
    };
    // Create database - train - align
    let lastdb = std::process::Command::new("lastdb")
        .args(&["-R", "00", "-Q", "0", &db_name, reference])
        .output()
        .ok()?;
    if !lastdb.status.success() {
        eprintln!("{:?} or {:?} might be invalid.", db_name, reference);
        panic!("lastdb-{}", String::from_utf8_lossy(&lastdb.stderr));
    }
    let p = format!("{}", thread);
    let last_train = std::process::Command::new("last-train")
        .args(&["-P", &p, "-Q", "0", &db_name, &reads])
        .output()
        .unwrap();
    if !last_train.status.success() {
        panic!("last-train-{}", String::from_utf8_lossy(&last_train.stderr));
    }
    let param = {
        let mut param = c_dir.clone();
        param.push("param.par");
        let mut wtr = BufWriter::new(std::fs::File::create(&param).unwrap());
        wtr.write_all(&last_train.stdout).unwrap();
        wtr.flush().unwrap();
        param.into_os_string().into_string().unwrap()
    };
    let lastal = std::process::Command::new("lastal")
        .args(&[
            "-f", "maf", "-P", &p, "-R", "00", "-Q", "0", "-p", &param, &db_name, &reads,
        ])
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("failed to invoke lastal");
    let raw_alignment = lastal.stdout.expect("failed to chatch lastal output");
    let last_split = std::process::Command::new("last-split")
        .stdin(std::process::Stdio::from(raw_alignment))
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("failed to invoke last-split");
    let raw_alignment = last_split.stdout.expect("failed to exec last-split");
    let maf_convert = std::process::Command::new("maf-convert")
        .args(&["tab", "--join", "1000"])
        .stdin(std::process::Stdio::from(raw_alignment))
        .output()
        .expect("failed to exec maf-convert");
    let alignments: Vec<_> = String::from_utf8_lossy(&maf_convert.stdout)
        .lines()
        .filter(|e| !e.starts_with("#"))
        .filter_map(|e| last_tiling::LastTAB::from_line(&e))
        .collect();
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir).ok()?;
    Some(alignments)
}

pub fn last_alignment<P: AsRef<std::path::Path>>(
    query: &P,
    refr: &P,
    threads: usize,
) -> Option<Vec<last_tiling::LastTAB>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 10_000;
    let mut c_dir = std::env::current_dir().ok()?;
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir).ok()?;
    // Create reference and reads.
    let (reference, reads) = (refr.as_ref().to_str()?, query.as_ref().to_str()?);
    let db_name = {
        let mut temp = c_dir.clone();
        temp.push("reference");
        temp.into_os_string().into_string().unwrap()
    };
    // Create database - train - align
    let lastdb = std::process::Command::new("lastdb")
        .args(&["-R", "00", "-Q", "0", &db_name, &reference])
        .output()
        .ok()?;
    if !lastdb.status.success() {
        panic!("lastdb-{}", String::from_utf8_lossy(&lastdb.stderr));
    }
    let p = format!("{}", threads);
    let lastal = std::process::Command::new("lastal")
        .args(&[
            "-f", "tab", "-P", &p, "-R", "00", "-Q", "0", &db_name, &reads,
        ])
        .output()
        .expect("failed to invoke lastal");
    let alignments: Vec<_> = String::from_utf8_lossy(&lastal.stdout)
        .lines()
        .filter(|e| !e.starts_with("#"))
        .filter_map(|e| last_tiling::LastTAB::from_line(&e))
        .collect();
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir).ok()?;
    Some(alignments)
}

pub fn dump_viewer(
    results: &HashMap<String, u8>,
    reads: &[EncodedRead],
    clusters: &[Cluster],
    contigs: &Contigs,
) -> std::io::Result<String> {
    let clusters = summarize_clusters(clusters, results);
    let summary = convert_to_d3_data_with_assign(&contigs, &reads, &clusters, &results);
    Ok(serde_json::ser::to_string(&summary).unwrap())
}

fn summarize_clusters(clusters: &[Cluster], results: &HashMap<String, u8>) -> Vec<Cluster> {
    let max = results.values().copied().max().unwrap_or(0) as usize;
    let mut aggregated: Vec<_> = (0..=max)
        .map(|id| Cluster {
            id,
            members: vec![],
        })
        .collect();
    use last_decompose::find_breakpoint::ReadClassify;
    // Put initial clusters into the most probable cluster.
    for cluster in clusters {
        // Collect the reads and clusters
        let reads: std::collections::HashMap<u8, usize> = results
            .iter()
            .filter(|&(ref id, _)| cluster.has(id))
            .map(|e| e.1)
            .fold(HashMap::new(), |mut x, &y| {
                *x.entry(y).or_default() += 1;
                x
            });
        // Determine which clusters these initial cluster should be put on.
        if let Some((&argmax, _)) = reads.iter().max_by_key(|x| x.1) {
            let members = cluster.members.iter().cloned().map(|mut mem| {
                mem.cluster = argmax as usize;
                mem
            });
            aggregated[argmax as usize].members.extend(members);
        }
    }
    aggregated
}
