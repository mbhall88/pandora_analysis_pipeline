import shutil
from pathlib import Path

rule add_denovo_paths_to_msa:
    input:
        denovo_dirs = expand("analysis/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery/denovo_paths", sample=config["samples"]),
        msa = "data/msas/{clustering_tool}/{gene}.fa"
    output:
        "analysis/{coverage}x/msas/{clustering_tool}/{gene}.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 200 * attempt
    #group: "add_denovo"
    run:
        new_msa_path = Path(output[0])
        old_msa_path = Path(input.msa)

        if not new_msa_path.parent.is_dir():
            new_msa_path.parent.mkdir(parents=True, exist_ok=True)

        denovo_paths = []
        for denovo_dir in input.denovo_dirs:
            denovo_paths.extend(list(Path(denovo_dir).glob(f"{wildcards.gene}.*.fa")))

        if not denovo_paths:
            shutil.copy(old_msa_path, new_msa_path)
        else:
            with new_msa_path.open("w") as fh_out:
                fh_out.write(old_msa_path.read_text())

                for p in denovo_paths:
                    read_counter = 1
                    sample = p.parts[-4]
                    name = p.with_suffix("").name

                    with p.open() as fasta:
                        for line in fasta:
                            if line.startswith(">"):
                                fh_out.write(f"{line.rstrip()}_sample={sample}_{name}_path{read_counter}\n")
                                read_counter += 1
                            else:
                                fh_out.write(line)

rule run_msa_after_adding_denovo_paths:
    input:
        updated_msa = "analysis/{coverage}x/msas/{clustering_tool}/{gene}.fa",
        original_msa = "data/msas/{clustering_tool}/{gene}.fa"
    output:
        "analysis/{coverage}x/msas/{clustering_tool}/{gene}.clustalo.fa"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000
    singularity:
        "docker://quay.io/biocontainers/clustalo:1.2.4--hfc679d8_2"
    #group: "add_denovo"
    log:
        "logs/{coverage}x/run_msa_after_adding_denovo_paths/{clustering_tool}/{gene}.log"
    shell:
        """
        if cmp -s {input.updated_msa} {input.original_msa}  # if files are the same
        then
            echo "Dont need to run MSA for this gene."
            cp {input.updated_msa} {output} > {log} 2>&1  # avoid having to rerun MSA
        else
            echo "Denovo paths have been added to this MSA. Running clustal omega..."
            clustalo --dealign --threads {threads} --in {input.updated_msa} --out {output} > {log} 2>&1
        fi
        """

rule build_prg_after_adding_denovo_paths:
    input:
        msa = "analysis/{coverage}x/msas/{clustering_tool}/{gene}.clustalo.fa",
        original_prg = "data/prgs/ecoli_pangenome_PRG_210619.fa",
        original_msa = "data/msas/{clustering_tool}/{gene}.fa", 
    output:
        prg = "analysis/{coverage}x/prgs/{clustering_tool}/{gene}.prg.fa",
    shadow: "shallow"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: (1000 + (attempt * 1000)) * attempt
        #mem_mb = lambda wildcards, attempt: 32000 * attempt
    params:
        script = "scripts/make_prg_from_msa.py",
        max_nesting_lvl = config.get("max_nesting_lvl", 5),
        prefix = lambda wildcards, output: output.prg.replace("".join(Path(output.prg).suffixes), "")
    #singularity: "containers/make_prg.simg"
    singularity: CONDA_IMG
    conda:
        "../envs/make_prg.yaml"
    #group: "add_denovo"
    log:
        "logs/{coverage}x/build_prg_after_adding_denovo_paths/{clustering_tool}/{gene}.log"
    shell:
        """
        if cmp -s {input.msa} {input.original_msa}  # if files are the same
        then  # dont run make_prg - get original prg entry for this gene
            echo "No denovo paths have been added. Getting original PRG."
            grep -A 1 {wildcards.gene} {input.original_prg} > {output.prg} 2> {log}
        else
            echo "Denovo paths have been added for this gene. Running make_prg"
            python3 {params.script} -v \
                --max_nesting {params.max_nesting_lvl} \
                --prefix {params.prefix} {input.msa} > {log} 2>&1
            mv {params.prefix}.max_nest{params.max_nesting_lvl}.min_match7.prg {output.prg} >> {log} 2>&1
            tmp_fname=$(mktemp)
            echo '>{wildcards.gene}' | awk 1 - {output.prg} > $tmp_fname && mv $tmp_fname {output.prg} >> {log} 2>&1
        fi
        """

prgs_file_names_after_denovo = []
for clustering_tool, gene in tool_msa_pair:
    for coverage in config["coverages"]:
        fname = f"analysis/{coverage}x/prgs/{clustering_tool}/{gene}.prg.fa"
        prgs_file_names_after_denovo.append(fname)

rule combine_prgs_after_adding_denovo_paths:
    input:
        prgs = prgs_file_names_after_denovo,
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
