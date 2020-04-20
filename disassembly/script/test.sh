#!/bin/bash
hai=TEST
test=$(cargo run --release --bin test -- ${hai})
echo ${test} "Is done"
