from pathlib import Path
import itertools

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


# ======================================================
# Functions
# ======================================================


# ======================================================
# Variables
# ======================================================
CONDA_IMG = config["conda_img"]
msa_paths = list(Path("data/msas").rglob("*.fa"))
cluster_name_to_path_map = {p.name.replace(".fa", ""): p for p in msa_paths}
assert len(msa_paths) == len(cluster_name_to_path_map)

tool_msa_pair = []
for gene_name in cluster_name_to_path_map:
    if gene_name.startswith("GC"):
        tool_msa_pair.append(("panx", gene_name))
    elif gene_name.startswith("Clus"):
        tool_msa_pair.append(("piggy", gene_name))

output_files = []

for sample, coverage, strategy in itertools.product(config["samples"], config["coverages"], config["subsample"]["strategies"]):
    output_files.extend([
        f"data/{sample}/{sample}.{coverage}x.{strategy}.nanopore.fastq"
        # f"analysis/{coverage}x/plots/{sample}/NanoPlot-report.html",
        # f"analysis/plots/{sample}/NanoPlot-report.html",
        # f"analysis/{coverage}x/compare_without_denovo/pandora_multisample_genotyped.vcf",
        # f"analysis/{coverage}x/{sample}/map_with_discovery/pandora_genotyped.vcf",
    ])

#
# output_files.extend(prgs_file_names_after_denovo)

# ======================================================
# Rules
# ======================================================
rule all:
    input:
        output_files

rules_dir = Path("rules/")
include: str(rules_dir / "filter.smk")
# include: str(rules_dir / "plot.smk")
# include: str(rules_dir / "index.smk")
# include: str(rules_dir / "map.smk")
# include: str(rules_dir / "create_new_prg.smk")
# include: str(rules_dir / "compare.smk")
