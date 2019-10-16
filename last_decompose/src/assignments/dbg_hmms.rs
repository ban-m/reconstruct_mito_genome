use last_tiling::Contigs;
use last_tiling::EncodedRead;
/// de Bruijn Graph powered by Hidden Markov model
/// It there is one markov model to each unit on contig.
#[derive(Debug, Default)]
pub struct DeBruijnGraphHiddenMarkovs<'a> {
    contigs: Vec<usize>,
    hmms: Vec<Vec<HiddenMarkov<'a>>>,
    param: Parameters,
    num_of_clusters: usize,
}

impl<'a> DeBruijnGraphHiddenMarkovs<'a> {
    // TODO
    pub fn new(_contigs: &Contigs, _num: usize) -> Self {
        Self::default()
    }
    // TODO
    pub fn push(&mut self, _color: usize, _read: &EncodedRead) {}
    // TODO
    pub fn predict(&self, _read: &EncodedRead) -> Vec<f64> {
        vec![]
    }
}

#[derive(Debug, Default)]
pub struct HiddenMarkov<'a> {
    seq: Vec<(usize, &'a [u8])>,
}

/// Global parameters.
#[derive(Debug, Default)]
pub struct Parameters {
    match_prob: f64,
    ins_start_prob: f64,
    del_start_prob: f64,
    ins_extend_prob: f64,
    del_extend_prob: f64,
    ins_to_del_prob: f64,
    del_to_ins_prob: f64,
}
