#!/bin/awk
BEGIN{
    OFS="\t"
}
($1 ~ /^@/){next} # Header file
(and($2,0x900) != 0x00) { next } # secondary alignment.
($3 ~ /NC_037304.1/){ print $1,"mitochondria"; next }
($3 ~ /NC_000932.1/){ print $1,"chroloplast"; next }
($3 ~ /00307/){print $1, "genome"; next}
{print $1,"*"}

