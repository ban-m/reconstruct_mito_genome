/// Definition of contigs.
use bio_utils::fasta;
use std::collections::HashMap;

/// A struct to represent contigs.
/// The value is [template, recvomp].
#[derive(Deserialize, Serialize, Debug, Default)]
pub struct Contigs {
    contigs: HashMap<String, [String; 2]>,
    names: Vec<String>,
}

impl std::fmt::Display for Contigs {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Names:")?;
        for (key, val) in &self.contigs {
            writeln!(f, ">{}( {} length)", key, val[0].as_bytes().len())?;
        }
        Ok(())
    }
}

impl Contigs {
    pub fn from_file<P: AsRef<std::path::Path>>(file: P) -> std::io::Result<Self> {
        fasta::parse_into_vec(file).map(Self::new)
    }
    pub fn new(records: Vec<fasta::Record>) -> Self {
        let contigs: HashMap<_, _> = records
            .iter()
            .map(|e| {
                let id = e.id().to_string();
                let template = String::from_utf8(e.seq().to_vec()).unwrap();
                let revcmp = String::from_utf8(revcmp(e.seq())).unwrap();
                (id, [template, revcmp])
            })
            .collect();
        let names: Vec<_> = records.iter().map(|e| e.id().to_string()).collect();
        Self { contigs, names }
    }
    pub fn get(&self, key: &str) -> Option<&[u8]> {
        self.contigs.get(key).map(|e| e[0].as_bytes())
    }
    pub fn get_revcmp(&self, key: &str) -> Option<&[u8]> {
        self.contigs.get(key).map(|e| e[1].as_bytes())
    }
    pub fn get_id(&self, key: &str) -> Option<u16> {
        self.names
            .iter()
            .enumerate()
            .filter_map(|(idx, name)| if name == key { Some(idx as u16) } else { None })
            .nth(0)
    }
    pub fn get_by_id(&self, id: u16) -> Option<&[u8]> {
        self.names.get(id as usize).and_then(|e| self.get(e))
    }
    pub fn get_by_id_revcmp(&self, id: u16) -> Option<&[u8]> {
        self.names.get(id as usize).and_then(|e| self.get_revcmp(e))
    }
    pub fn names(&self) -> &[String] {
        &self.names
    }
}

#[inline]
fn revcmp(seq: &[u8]) -> Vec<u8> {
    seq.into_iter()
        .rev()
        .map(|&e| match e {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => unreachable!(),
        })
        .collect()
}
