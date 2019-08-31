# cargo run --release --bin main -- \
#       /grid/ban-m/arabidopsis_thaliana/sequel/dilution_onthefly2/selected_2nd.bam \
#       /grid/ban-m/arabidopsis_thaliana/sequel/dilution_onthefly2/selected_2nd.fq

      
cargo run --release --bin alignment_recover -- \
      ./temp/ovlp_cntn_trns_hvst.sam \
      /grid/ban-m/arabidopsis_thaliana/sequel/dilution_onthefly2/selected_2nd.fq
