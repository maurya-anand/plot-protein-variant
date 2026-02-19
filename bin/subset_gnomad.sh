#!/usr/bin/env bash

TRANSCRIPT_TSV=$1

GNOMAD_CSV=$2

seq_region_name=$(awk -F'\t' '$1 == "seq_region_name" {print $2}' ${TRANSCRIPT_TSV})

start=$(awk -F'\t' '$1 == "start" {print $2}' ${TRANSCRIPT_TSV})

end=$(awk -F'\t' '$1 == "end" {print $2}' ${TRANSCRIPT_TSV})

LC_ALL=C awk -F',' -v chrom="$seq_region_name" -v start="$start" -v end="$end" '
    BEGIN {OFS=","}
    NR == 1 || ($1 == chrom && $2 >= start && $2 <= end)
' ${GNOMAD_CSV} > subset.csv

realpath subset.csv