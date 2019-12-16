#!/usr/bin/env bash
set -eux

pandora_executable="$1"
prg_file=$(realpath "$2")
input_reads=$(realpath "$3")
outdir=$(realpath "$4")
threads="$5"
technology="$6"
log_level="$7"
log=$(realpath "$8")
use_discover="$9"
input_ref="$10"

if [ "${use_discover,,}" = "true" ]; then
    discover="--discover"
elif [ "${use_discover,,}" = "false" ]; then
    discover=" "
else
    echo "Error: Unrecognised option passed for use_discover. Must be true or false. Got: ${use_discover}"
    exit 1
fi

if [ "${technology,,}" = "illumina" ]; then
    technology_param="--illumina"
elif [ "${technology,,}" = "nanopore" ]; then
    technology_param=" "
else
    echo "Error: Unrecognised option passed for technology. Must be illumina or nanopore. Got: ${technology}"
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
    --max_covg 100000 \
    ${technology_param} \
    "$discover" > "$log" 2>&1
