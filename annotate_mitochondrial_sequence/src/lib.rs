#![allow(non_snake_case)]
#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
// Species Name -> (start, end, genename, strand)
// where strand is -1 if reverse, 1 if foward.
#[derive(Debug)]
pub struct Annotation {
    pub gene: String,
    pub start: usize,
    pub end: usize,
    pub strand: i8,
}

impl Annotation {
    fn new(input: &str) -> Option<(String, Self)> {
        let contents: Vec<&str> = input.split('\t').collect();
        let spname = contents[0].to_string();
        let gene = contents[1].to_string();
        let start: usize = contents[2].parse().ok()?;
        let end: usize = contents[3].parse().ok()?;
        let strand = if contents[4] == "+" { 1 } else { -1 };
        let res = Annotation {
            gene,
            start,
            end,
            strand,
        };
        Some((spname, res))
    }
}

pub fn open_annotation(file: &str) -> std::io::Result<HashMap<String, Vec<Annotation>>> {
    let mut res = HashMap::new();
    for (spname, annot) in BufReader::new(File::open(&Path::new(file))?)
        .lines()
        .skip(1)
        .filter_map(|e| e.ok())
        .filter_map(|e| Annotation::new(&e))
    {
        let entry = res.entry(spname).or_insert_with(Vec::new);
        entry.push(annot);
    }
    for (_spname, ants) in res.iter_mut() {
        ants.sort_by_key(|e| e.start)
    }
    Ok(res)
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Map {
    contig_name: String,
    is_tRNA: bool,
    is_forward: bool,
    start: u64,
    end: u64,
    ref_species: String,
    gene_name: String,
    confidence: f64,
}

impl Map {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        contig_name: &str,
        is_tRNA: bool,
        is_forward: bool,
        start: u64,
        end: u64,
        ref_species: &str,
        gene_name: &str,
        confidence: f64,
    ) -> Self {
        let contig_name = contig_name.to_string();
        let ref_species = ref_species.to_string();
        let gene_name = gene_name.to_string();
        Map {
            contig_name,
            is_tRNA,
            is_forward,
            start,
            end,
            ref_species,
            gene_name,
            confidence,
        }
    }
}
