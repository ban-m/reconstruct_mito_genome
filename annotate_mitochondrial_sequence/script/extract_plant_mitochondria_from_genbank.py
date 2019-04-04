"""
A Modulo to convert genbank file into fasta and tab-delimited feature files.
Input: python3.6 [this script] [GenBank file] [file to write fasta] [file to write feature]
"""
import sys
import itertools as it
import Bio
from Bio import SeqIO

def is_mitochondrial(record):
    """
    determine whether given records is derived from
    plant mitochondria. 
    In this case, I define 'plant' as Angiosperms.
    It is equivalent to use 'Magnoliophyta' as key.
    Input: a genbank record
    Output: bool
    """
    return 'Magnoliophyta' in record.annotations['taxonomy']

def extract_fasta(record):
    """
    Extract fasta sequence from given record.
    Input: a genbank record.
    Output: fasta record
    """
    seqid = record.annotations['organism'].replace(' ', '_')
    seq = record.seq
    description = "converted by BM from GenBank at 20190403"
    rec = Bio.SeqRecord.SeqRecord(seq, seqid, description=description)
    return rec

def feature_to_tuple(feature):
    """
    Convert a feature to (gene name, start, end, strand)
    Input: a SeqFeature
    Ouput: (str, Int, Int)
    """
    start = feature.location.parts[0].start.position + 1
    end = feature.location.parts[-1].end.position
    # Maybe the name is not available.
    if 'gene' in feature.qualifiers:
        name = feature.qualifiers['gene'][0]
    else:
        name = next(filter(lambda x: x.startswith("GeneID"), feature.qualifiers['db_xref']))
    strand = feature.strand
    if strand == 1:
        return (name, start, end, '+')
    return (name, start, end, '-')

def extract_feature(record):
    """
    Extract features as tuple of (species name, gene name, start, end, strand).
    This is sufficient because the genome has only one contig, usually.
    Input: a genbank record.
    Output: a iterator over (String, String, Int, Int, char)
    """
    spname = record.annotations['organism'].replace(' ', '_')
    return map(lambda tup: (spname,) + tup,
               map(feature_to_tuple,
                   filter(lambda x: x.type == 'gene',
                          record.features)))

def filter_plant_mitochondria(genbank, fasta_wtr, feature_wtr):
    """
    Picking plant mitochondrial records from given input, then
    Output fasta file and feature file.
    Input: genbank reader(generator), file handler to fasta records,
    and file handler to feature records.
    """
    # Should be a list. Not iterable.
    genbank = list(filter(is_mitochondrial, genbank))
    # genbank = list(genbank)
    # genbank would not be consumed
    seqs = map(extract_fasta, genbank)
    Bio.SeqIO.write(seqs, fasta_wtr, 'fasta-2line')
    print("SpeciesName\tGeneNameOrID\tStart\tEnd\tStrand", file=feature_wtr)
    feature_tuples = it.chain.from_iterable(map(extract_feature, genbank))
    for (species, genename, start, end, strand) in feature_tuples:
        print("{}\t{}\t{}\t{}\t{}"
              .format(species, genename, start, end, strand),
              file=feature_wtr)

def main(args):
    """
    Main Function to parse GenBank format and output plant mitochondrial DNA(in fasta format and tab-delimited file).
    """
    with open(args[0], 'r') as genbank:
        with open(args[1], 'w') as fasta_handle:
            with open(args[2], 'w') as feature_wtr:
                genbank = Bio.SeqIO.parse(genbank, 'genbank')
                filter_plant_mitochondria(genbank, fasta_handle, feature_wtr)


if __name__ == '__main__':
    ARGS = sys.argv;
    if len(ARGS) != 4:
        print('Error: Please supply exactly three argument.', file=sys.stderr)
        exit(1)
    main(ARGS[1:4])
