use bio_utils::fasta::Record;
use last_tiling::LastTAB;
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Summary {
    references: Vec<Reference>,
    alignments: Vec<Alignment>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Reference {
    id: String,
    length: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Alignment {
    reference_name: String,
    ref_start: usize,
    ref_end: usize,
    query_name: String,
    query_start: usize,
    query_end: usize,
}
pub fn annotate_aln_contigs_to_ref(refs: &[Record], alns: &[LastTAB], thr: usize) -> Summary {
    let alignments: Vec<_> = alns
        .into_iter()
        .filter(|aln| (aln.seq1_matchlen() + aln.seq2_matchlen()) / 2 > thr)
        .map(|aln| {
            let reference_name = aln.seq1_name().to_string();
            let ref_start = aln.seq1_start_from_forward();
            let ref_end = aln.seq1_end_from_forward();
            let query_name = aln.seq2_name().to_string();
            let query_start = aln.seq2_start_from_forward();
            let query_end = aln.seq2_end_from_forward();
            if aln.seq1_direction() == aln.seq2_direction() {
                Alignment {
                    reference_name,
                    ref_start,
                    ref_end,
                    query_name,
                    query_start,
                    query_end,
                }
            } else {
                Alignment {
                    reference_name,
                    ref_start,
                    ref_end,
                    query_name,
                    query_start: query_end,
                    query_end: query_start,
                }
            }
        })
        .collect();
    let references: Vec<_> = refs
        .into_iter()
        .map(|rec| {
            let id = rec.id().to_string();
            let length = rec.seq().len();
            Reference { id, length }
        })
        .collect();
    Summary {
        references,
        alignments,
    }
}
