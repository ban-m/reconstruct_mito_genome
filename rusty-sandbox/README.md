# Rusty sandbox

This is a experimental tools under development.

# Tools

## alignment_recover.rs

This binary reconstruct an alignment from specified sam file and reference fastq records.
Currently, this binary is for all vs all alignment. Thus, the input format is as follows:
```bash
cargo run --release --bin alignment_recover -- [sam] [fastq]
```

## filter_contained_read.rs

This binary is to clean a sam file. This filter trims all the contained alignment. In other
words, it detects all alignment which completely consumes either query or reference and remove them from
sam file. It outputs also sam file.

```bash
cargo run --release --bin filter_contained_read -- [sam] [fastq] [sam-output] [fastq-output]
```

## decompose_into_connected_component.rs

This program is to calculate how many connected components are there in sam file.
```bash
cargo run --release --bin decompose_into_connected_components -- [sam file]
```
The output format is tab-delimited text file. The fields are [component number], [# of node in the component], and [# of edges in the component]. Note that the number begins from 2. This is rather a programming convention, not a bug.


## transitive_edge_reduction.rs

This program removes all transitive edges from input sam file.
Inside this program, it first converts the sam file into a directed graph, then
detects all edges which is transitive.
```bash
cargo run --release --bin transitive_edge_reduction -- [sam]
```

## decompose_into_connected_st_component.rs

This program decomposes the input sam file into strongly connected components, and output the stats.
```bash
cargo run --release --bin decompose_into_st_components -- [sam]
```
The output format is the same as that of decompose_into_connected_componet, the undirected version.

## select_heviest_edge.rs

This program selectes an edges for each node out of all the out edges from that node, and remove the rest edges.
This reduction is too strong to be used.
```bash
cargo run --release --bin select_heviest_edge [sam]
```

## error_rate_of_pileup.rs

This program compute the error profile from all vs all sam file. Currently ,the sample rate is 0.001.
```bash
cargo run --release --bin error_rate_of_fileup -- [sam] [fastq]
```
The output format is tsv, [position], [substitution number], [insertion bases], [deletion bases], and [sequencing depth]. Note that the depth is defined as the sum of match. So, the insertion and deletion can be exceed the depth. Insertion is the insertion after the reference base.


## select_overlaps.rs
This program retains all the overlaps, while triming all the `partial` alignments.
Specifically, it require  `(ref_begin < THR && query_end < THR) || (query_begin < THR && ref_end < THR)`,
where `ref` and `query` means reference position and query position, respectively. So, if I use

...... as proper alignments, and """""" as clipped region, the following two are accepted.


Reference:    ...............""""""""""
Query    :""""...............

Or

Reference:""""""""...................
Query    :        ..................."""""""""

```bash
cargo run --release --bin select_overlaps --- [sam]
```

## output_pileup.rs

This is a  captive program to calculate the pileup and output to a specified file.
```bash
cargo run --release --bin output_puleup -- [sam] [fastq] 
```

Then,
```bash
$> 1243 # Name of a read you want to calculate the pileup.
$> ./temp # filename to which the output would be written.
```
