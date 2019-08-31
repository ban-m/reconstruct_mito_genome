# Tiling analysis




## Input format

File := Empty | Record + "\n"
Empty := ""
i.e., The input file is a line-by-line format.

Record = [Read ID] + ":" + Encoded
Encoded = Empty | " " + Unit + Encoded
Unit = "G-" + Int | Int + "-" + Int

