#!/usr/bin/env bash
set -eux

pandora_executable="$1"
prg_file=$(realpath "$2")
input_reads=$(realpath "$3")
outdir=$(realpath "$4")
threads="$5"
log_level="$6"
log=$(realpath "$7")
use_discover="$8"
input_ref="$9"

if [ "${use_discover,,}" = "true" ]; then
    discover="--discover"
elif [ "${use_discover,,}" = "false" ]; then
    discover=" "
else
    echo "Error: Unrecognised option passed for use_discover. Must be true or false. Got: ${use_discover}"
    exit 1
fi

genome_size=$(grep -v '>' "$input_ref" | wc | awk '{ print $3-$1 }')
mkdir -p "$outdir"
cd "$outdir" || exit 1

"$pandora_executable" map --prg_file "$prg_file" \
    --read_file "$input_reads" \
    --outdir "$outdir" \
    -t "$threads" \
    --output_kg \
    --genome_size "$genome_size" \
    --output_covgs \
    --output_vcf \
    --log_level "$log_level" \
    --genotype \
    "$discover" > "$log" 2>&1
