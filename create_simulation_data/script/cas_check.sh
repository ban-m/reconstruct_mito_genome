#!/bin/bash

function check() {
    echo $1
    tail -n+2 $1 | awk '{sum +=1; ok += $2 ==$3 }END{print sum,ok, ok/sum}'
}


for file in ./result/cas*
do
    check $file
done
