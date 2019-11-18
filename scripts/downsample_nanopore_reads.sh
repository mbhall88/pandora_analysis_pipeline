#!/usr/bin/env bash
if [[ "$#" -lt 5 ]]; then
    echo "Error: Illegal number of parameters"
    echo -e "usage:\n$(basename "$0") <reads> <ref> <covg> <outname> <strategy> [<min_len> <q_weight>]"
    exit 1
fi

set -euv

reads="$1"
ref="$2"
covg="$3"
outname="$4"
strategy="$5"
min_len="${6:-1}"
q_weight="${7:-10}"

genome_size=$(grep -v '^>' "$ref" | wc | awk '{print $3-$1}')
num_bases_to_keep=$((genome_size * covg))

if [ "$strategy" = "filter" ]
then
    filtlong \
        --verbose \
        --target_bases "$num_bases_to_keep" \
        --min_length "$min_len" \
        --mean_q_weight "$q_weight" \
        "$reads" > "$outname"
elif [ "$strategy" = "random" ]
then
    rasusa --input "$reads" \
        --coverage "$covg" \
        --genome-size "$genome_size" \
        --output "$outname" \
        -v
else
    echo "Invalid strategy given: $strategy"
    exit 1
fi
