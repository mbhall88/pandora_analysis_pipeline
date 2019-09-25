from pathlib import Path


rule add_denovo_paths:
    input:
        denovo_dirs = expand("analysis/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery/denovo_paths", sample=config["samples"]),
        msa = "data/msas/{clustering_tool}/{gene}.fa"
    output:
        updated_msa = "analysis/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.clustalo.fa",
        appended_msa = "analysis/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.fa",
        prg = "analysis/{coverage}x/{sub_strategy}/updated_prgs/{clustering_tool}/{gene}.prg.fa"
    threads: 2
    shadow: "shallow"
    resources:
        mem_mb = lambda wildcards, attempt: {1: 1000, 2: 4000, 3: 16000}.get(attempt, 64000)
    params:
        log_level = "DEBUG",
        make_prg_script = "scripts/make_prg_from_msa.py",
        max_nesting_lvl = config.get("max_nesting_lvl", 5),
        prefix = lambda wildcards, output: output.prg.replace("".join(Path(output.prg).suffixes), "")
    singularity: CONDA_IMG
    conda:
        "../envs/add_denovo_paths.yaml"
    log:
        "logs/add_denovo_paths/{coverage}x/{sub_strategy}/{clustering_tool}/{gene}.log"
    script:
        "scripts/add_denovo_paths.py"


rule extract_original_prg:
    input:
        original_prg = config["original_prg"]
    output:
        prg = "analysis/{coverage}x/{sub_strategy}/original_prgs/{clustering_tool}/{gene}.prg.fa"
    threads: 1
    resources:
        mem_mb = 200
    log:
       "logs/extract_original_prg/{coverage}x/{sub_strategy}/{clustering_tool}/{gene}.log"
    shell:
        "grep -A 1 {wildcards.gene} {input.original_prg} > {output.prg} 2> {log}"


def aggregate_prgs_input(wildcards):
    checkpoint_output = checkpoints.map_with_discovery.get(coverage=wildcards.coverage,
                                                           sub_strategy=wildcards.sub_strategy
                                                           ).output.denovo_dir
    gene_tool_pairs_in_denovo = get_gene_tool_pairs_in_denovo_dir(checkpoint_output)
    input_files = []
    for tool, gene in TOOL_MSA_PAIR:
        if (gene, tool) in gene_tool_pairs_in_denovo:
            # get filename for add denovo paths
            input_files.append(f"analysis/{wildcards.coverage}x/{wildcards.sub_strategy}/updated_prgs/{tool}/{gene}.prg.fa")
        else:
            # get filename for extract prgs_in_original prg
            input_files.append(f"analysis/{wildcards.coverage}x/{wildcards.sub_strategy}/original_prgs/{tool}/{gene}.prg.fa")

    return input_files

def get_gene_tool_pairs_in_denovo_dir(denovo_dir: str) -> set:
    gene_tool_pairs = set()
    denovo_files = Path(denovo_dir).rglob("*discovery.fa")
    for path in denovo_files:
        name = path.name
        gene = name.split(".")[0]
        if gene.startswith("Cluster"):
            gene_tool_pairs.add((gene, "piggy"))
        elif gene.startswith("GC"):
            gene_tool_pairs.add((gene, "panx"))
        else:
            raise ValueError(f"Cannot find clustering tool for {name}")

    return gene_tool_pairs

rule aggregate_prgs:
    input:
        prgs = aggregate_prgs_input,
        original_prg = "data/prgs/ecoli_pangenome_PRG_210619.fa"
    output:
        "analysis/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 500 * attempt
    run:
        import fileinput
        with open(output[0], "w") as fout, fileinput.input(input.prgs) as fin:
            for line in fin:
                # strip and add newline in case some lines are missing newline
                fout.write(line.rstrip() + "\n")

        # check original prg and new prg have the same number of sequences
        prgs_in_original = 0
        with open(input.original_prg) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_original += 1

        prgs_in_new = 0
        with open(output[0]) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_new += 1

        assert prgs_in_original == prgs_in_new, "Original PRG and new PRG dont have the same number of entries!"
