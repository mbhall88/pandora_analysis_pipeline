#!/usr/bin/env bash
if [[ "$#" -lt 6 ]]; then
    echo "Error: Illegal number of parameters"
    echo -e "usage:\n$(basename "$0") <reads_1> <reads_2> <ref> <covg> <outname_1> <outname_2>"
    exit 1
fi

set -euv

reads_1="$1"
reads_2="$2"
ref="$3"
covg="$4"
outname_1="$5"
outname_2="$6"

genome_size=$(grep -v '^>' "$ref" | wc | awk '{print $3-$1}')
num_bases_to_keep=$((genome_size * covg))

python downsample_illumina_pe_reads/downsample_illumina_pe_reads.py --reads1 $reads_1 --reads2 $reads_2 \
--number_of_bases $num_bases_to_keep --out_reads1 $outname_1 --out_reads2 $outname_2