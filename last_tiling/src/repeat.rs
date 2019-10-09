//! The definitions of repeats.
use serde_json;

/// A repeat pairs.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RepeatPairs {
    reps: Vec<Repeat>,
}

/// A repeat. Note that the id should be consistent with other data such as contigs.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Repeat {
    id: u16,
    name: String,
    start: usize,
    end: usize,
}

/// Recover from path
pub fn open<P: AsRef<std::path::Path>>(file: P) -> std::io::Result<Vec<RepeatPairs>> {
    use std::io::BufReader;
    Ok(serde_json::de::from_reader(std::fs::File::open(file).map(BufReader::new)?).unwrap())
}

impl RepeatPairs {
    pub fn new(aln: &super::LastTAB, contigs: &super::Contigs) -> Option<Self> {
        let mut reps = Vec::with_capacity(2);
        reps.push(Repeat {
            id: contigs.get_id(aln.seq1_name())?,
            name: aln.seq1_name().to_string(),
            start: aln.seq1_start_from_forward(),
            end: aln.seq1_end_from_forward(),
        });
        reps.push(Repeat {
            id: contigs.get_id(aln.seq2_name())?,
            name: aln.seq2_name().to_string(),
            start: aln.seq2_start_from_forward(),
            end: aln.seq2_end_from_forward(),
        });
        Some(Self { reps })
    }
    pub fn inner(&self)->&[Repeat]{
        &self.reps
    }
    pub fn len(&self) -> usize {
        self.reps.len()
    }
    pub fn is_empty(&self) -> bool {
        self.reps.is_empty()
    }
}

impl Repeat {
    pub fn id(&self) -> u16 {
        self.id
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn start(&self) -> usize {
        self.start
    }
    pub fn end(&self) -> usize {
        self.end
    }
}

impl std::ops::Index<usize> for RepeatPairs {
    type Output = Repeat;
    fn index(&self, index: usize) -> &Self::Output {
        self.reps.index(index)
    }
}
