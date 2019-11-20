#!/bin/bash
function check() {
    echo $1
    tail -n+2 $1 | awk '{sum +=1; ok += $2 ==$3 }END{print sum,ok, ok/sum}'
}

check ./result/mock_genome/predictions_two_backs_with_merge.tsv
check ./result/mock_genome/predictions_two_backs_no_merge.tsv
check ./result/mock_genome/predictions_one_back_with_merge.tsv
check ./result/mock_genome/predictions_one_back_no_merge.tsv
