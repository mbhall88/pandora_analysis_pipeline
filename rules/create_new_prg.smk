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

prgs_file_names_after_denovo = []
for clustering_tool, gene in tool_msa_pair:
    for coverage in config["coverages"]:
        fname = f"analysis/{coverage}x/prgs/{clustering_tool}/{gene}.prg.fa"
        prgs_file_names_after_denovo.append(fname)

def aggregate_prgs_input(wildcards):
    checkpoint_output = checkpoints.map_with_discovery.get(coverage=wildcards.coverage,
                                                           sub_strategy=wildcards.sub_stractegy,
                                                           sample=wildcards.sample).output.denovo_dir
    gene_tool_pairs_in_denovo = get_gene_tool_pairs_in_denovo_dir(checkpoint_output)
    input_files = []
    for gene, tool in tool_msa_pair:
        if (gene, tool) in gene_tool_pairs_in_denovo:
            # get filename for add denovo paths
        else:
            # get filename for extract original prg

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
        "analysis/{coverage}x/prgs/denovo_updated.prg.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 500 * attempt
    log:
        "logs/{coverage}x/combine_prgs_after_adding_denovo_paths.log"
    run:
        for p in input.prgs:
            shell(f"awk 1 {p} >> {{output}} 2>> {{log}}")
        # check original prg and new combined one have the same number of sequences
        original = 0
        with open(input.original_prg) as fh:
            for line in fh:
                if line.startswith(">"):
                    original += 1

        combined = 0
        with open(output[0]) as fh:
            for line in fh:
                if line.startswith(">"):
                    combined += 1
        assert original == combined, "Original PRG and new combined PRG dont have the same number of entries!"
    #shell:
    #    """
    #    for prg in {input}
    #    do
    #        awk 1 $prg >> {output} 2> {log}
    #    done
    #    """
