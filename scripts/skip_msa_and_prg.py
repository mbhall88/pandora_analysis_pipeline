"""
This script is intended to reduce the number of jobs needed to be run to add denovo
candidate paths back into the PRG. It first determines which genes have candidate
paths after running denovo discovery. Then for any gene which does not have a candidate
path, it pulls out the original PRG entry for that gene from the original PRG and
writes it to file.
"""
from pathlib import Path
import argparse
import sys
import logging
import pysam
from typing import Dict, Set

LOGGING_LEVELS = {
    0: "NOTSET",
    1: "CRITICAL",
    2: "ERROR",
    3: "WARNING",
    4: "INFO",
    5: "DEBUG",
}
READ_CHUNKS = 250


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prg", required=True, help="""Path to the original PRG file.""", type=Path
    )
    parser.add_argument(
        "--coverage",
        required=True,
        help="""Search for denovo paths for this coverage.""",
        type=int,
    )
    parser.add_argument(
        "--log_level",
        help="Level of logging. 0 is none, 5 is for debugging. Default is 4 "
             "which will report info, warnings, errors, and critical information.",
        default=4,
        type=int,
        choices=range(6),
    )

    args = parser.parse_args()

    log_level = LOGGING_LEVELS.get(args.log_level)
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s]:%(levelname)s:%(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p ",
    )

    if not args.prg.is_file():
        raise FileNotFoundError(f"PRG file {args.prg} does not exist!")

    return args


def generate_denovo_search_pattern(coverage: int) -> str:
    return (
        f"analysis/{coverage}x/*/map_with_discovery/denovo_paths/*denovo_discovery.fa"
    )


def get_set_of_denovo_path_gene_names(coverage: int) -> Set[str]:
    denovo_paths = Path().rglob(generate_denovo_search_pattern(coverage))
    denovo_genes = set()
    for p in denovo_paths:
        gene = p.name.split(".")[0]
        denovo_genes.add(gene)

    return denovo_genes


def load_prg(prg_path: Path) -> Dict[str, str]:
    prg = dict()
    with pysam.FastxFile(prg_path) as fh:
        for record in fh:
            if record.name in prg:
                logging.warning(f"{record.name} is already in PRG dict.")
                continue
            prg[record.name] = record.sequence
    return prg


def generate_output_path_for_gene(gene: str, coverage: int) -> Path:
    if gene.startswith("GC"):
        clustering_tool = "panx"
    elif gene.startswith("Cluster"):
        clustering_tool = "piggy"
    else:
        raise ValueError(f"{gene} does not start with a recognised pattern.")

    return Path(f"analysis/{coverage}x/prgs/{clustering_tool}/{gene}.prg.fa")


def main():
    args = cli()

    logging.info("Getting denovo path gene names.")
    denovo_genes = get_set_of_denovo_path_gene_names(args.coverage)
    logging.info(f"Got {len(denovo_genes)} denovo genes.")

    logging.info("Loading PRG file.")
    prg = load_prg(args.prg)
    logging.info("PRG file loaded.")

    # go through all genes, if not in denovo path genes create prg file
    genes_not_in_denovo = prg.keys() - denovo_genes
    logging.info(f"{len(genes_not_in_denovo)} genes to process...")

    counter = 0
    for i, gene in enumerate(genes_not_in_denovo):
        logging.debug(f"Processing {gene}")
        sequence = prg[gene]
        gene_outpath = generate_output_path_for_gene(gene, args.coverage)

        if not gene_outpath.parent.is_dir():
            gene_outpath.parent.mkdir(parents=True, exist_ok=True)

        with gene_outpath.open("w") as fout:
            fout.write(f">{gene}\n{sequence}\n")

        counter += 1
        if counter == READ_CHUNKS:
            logging.info(f"Processed {i + 1} genes.")
            counter = 0


if __name__ == "__main__":
    main()