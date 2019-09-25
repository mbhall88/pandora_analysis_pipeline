import logging

from pathlib import Path
from snakemake.shell import shell
from typing import List, TextIO

log_level = snakemake.params.log_level
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)


def get_denovo_path_filepaths(denovo_dirs: List[str]) -> List[Path]:
    logging.debug("Getting filepaths for all de novo paths.")
    denovo_paths = []
    for denovo_dir in denovo_dirs:
        denovo_paths.extend(
            list(Path(denovo_dir).glob(f"{snakemake.wildcards.gene}.*.fa"))
        )

    logging.debug(f"Got {len(denovo_paths)} de novo paths")
    return denovo_paths


def append_denovo_paths_to_msa(
    denovo_paths: List[Path], fh_out: TextIO, old_msa: Path
) -> None:
    logging.info("Appending de novo paths to old MSA.")
    fh_out.write(old_msa.read_text())

    for p in denovo_paths:
        read_counter = 1
        sample = p.parts[-4]
        name = p.with_suffix("").name

        with p.open() as fasta:
            for line in fasta:
                if line.startswith(">"):
                    fh_out.write(
                        f"{line.rstrip()}_sample={sample}_{name}_path{read_counter}\n"
                    )
                    read_counter += 1
                else:
                    fh_out.write(line)

    logging.info("De novo paths appended to old MSA.")


def run_msa_after_adding_denovo_paths(
    appended_msa: str, updated_msa: str, threads: int
):
    logging.info("Running multiple sequence alignment.")
    shell(
        f"clustalo --dealign --threads {threads} --in {appended_msa} --out {updated_msa}"
    )
    logging.info("Multiple sequence alignment finished.")


def build_prg_after_adding_denovo_paths(
    make_prg_script: str,
    max_nesting_lvl: int,
    prefix: str,
    msa: str,
    prg: TextIO,
    gene: str,
):
    logging.info("Building PRG for MSA.")
    shell(
        f"python3 {make_prg_script} -v --max_nesting {max_nesting_lvl} --prefix {prefix} {msa}"
    )
    logging.info("Finished building PRG.")
    logging.debug("Adding header info to PRG and renaming to correct filepath.")
    prg.write(f">{gene}\n")
    tmp_prg = Path(f"{prefix}.max_nest{max_nesting_lvl}.min_match7.prg")
    prg.write(tmp_prg.read_text() + "\n")


def main():
    appended_msa = Path(snakemake.output.appended_msa)
    updated_msa = Path(snakemake.output.updated_msa)
    old_msa = Path(snakemake.input.msa)
    prg = Path(snakemake.output.prg)
    denovo_dirs = snakemake.input.denovo_dirs

    if not appended_msa.parent.is_dir():
        logging.debug("Creating parent directory for new MSA.")
        appended_msa.parent.mkdir(parents=True, exist_ok=True)

    denovo_paths = get_denovo_path_filepaths(denovo_dirs)

    if not denovo_paths:
        raise Exception(f"No de novo paths found for {snakemake.wildcards.gene}")

    with appended_msa.open("w") as fh_out:
        append_denovo_paths_to_msa(denovo_paths, fh_out, old_msa)

    run_msa_after_adding_denovo_paths(
        str(appended_msa), str(updated_msa), snakemake.threads
    )
    with prg.open("w") as prg_fh:
        build_prg_after_adding_denovo_paths(
            snakemake.params.make_prg_script,
            snakemake.params.max_nesting_lvl,
            snakemake.params.prefix,
            str(updated_msa),
            prg_fh,
            snakemake.wildcards.gene,
        )


main()
