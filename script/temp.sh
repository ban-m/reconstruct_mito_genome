#!/bin/bash
set -ue
RUST_LOG=debug cargo run --release --bin test_last_decompose_multiple -- \
        150 0 123 0.001 0.4 0.5 2> ./logfiles/test.log
