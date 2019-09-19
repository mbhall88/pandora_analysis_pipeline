#!/usr/bin/env bash
if [[ "$#" -lt 4 ]]; then
    echo "Error: Illegal number of parameters"
    echo -e "usage:\n$(basename "$0") <reads> <ref> <covg> <outname> [<min_len> <q_weight> <img>]"
    exit 1
fi

set -euv

reads="$1"
ref="$2"
covg="$3"
outname="$4"
min_len="${5:-1}"
q_weight="${6:-10}"

genome_size=$(grep -v '^>' "$ref" | wc | awk '{print $3-$1}')
num_bases_to_keep=$((genome_size * covg))

filtlong \
    --verbose \
    --target_bases "$num_bases_to_keep" \
    --min_length "$min_len" \
    --mean_q_weight "$q_weight" \
    "$reads" > "$outname"
