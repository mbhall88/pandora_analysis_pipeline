from pathlib import Path
import itertools
from snakemake.utils import min_version, validate

min_version("5.4.0")

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

TOOL_MSA_PAIR = []
for gene_name in cluster_name_to_path_map:
    if gene_name.startswith("GC"):
        TOOL_MSA_PAIR.append(("panx", gene_name))
    elif gene_name.startswith("Clus"):
        TOOL_MSA_PAIR.append(("piggy", gene_name))

output_files = []

for sample, coverage, strategy in itertools.product(config["samples"],
                                                    config["coverages"],
                                                    config["subsample"]["strategies"]):
    output_files.extend([
        f"analysis/{coverage}x/{strategy}/{sample}/map_with_discovery/pandora_genotyped.vcf",
        f"analysis/{coverage}x/{strategy}/prgs/denovo_updated.prg.fa",
        f"analysis/{coverage}x/{strategy}/compare_no_denovo/pandora_multisample_genotyped.vcf",
        f"analysis/{coverage}x/{strategy}/compare_with_denovo/pandora_multisample_genotyped.vcf",
        # f"analysis/{coverage}x/plots/{sample}/NanoPlot-report.html",
        # f"analysis/plots/{sample}/NanoPlot-report.html",
    ])

# ======================================================
# Rules
# ======================================================
rule all:
    input:
         output_files

rules_dir = Path("rules/")
include: str(rules_dir / "subsample.smk")
# include: str(rules_dir / "plot.smk")
include: str(rules_dir / "index.smk")
include: str(rules_dir / "map.smk")
include: str(rules_dir / "create_new_prg.smk")
include: str(rules_dir / "compare.smk")
