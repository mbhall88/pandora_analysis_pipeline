from pathlib import Path


rule add_denovo_paths:
    input:
        map_with_discovery_dirs = expand("analysis/{{technology}}/{{coverage}}x/{{sub_strategy}}/{sample}/map_with_discovery", sample=config["samples"]),
        msa = "data/msas/{clustering_tool}/{gene}.fa"
    output:
        updated_msa = "analysis/{technology}/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.clustalo.fa",
        appended_msa = "analysis/{technology}/{coverage}x/{sub_strategy}/msas/{clustering_tool}/{gene}.fa",
        prg = "analysis/{technology}/{coverage}x/{sub_strategy}/prgs/{clustering_tool}/{gene}.prg.fa"
    threads: 2
    shadow: "shallow"
    resources:
        mem_mb = lambda wildcards, attempt: {1: 1000, 2: 4000, 3: 16000}.get(attempt, 64000)
    params:
        log_level = "DEBUG",
        make_prg_script = "scripts/make_prg_from_msa.py",
        max_nesting_lvl = config.get("max_nesting_lvl", 5),
        prefix = lambda wildcards, output: output.prg.replace("".join(Path(output.prg).suffixes), ""),
        original_prg = config["original_prg"],
        denovo_dirs = lambda wildcards, input: [map_with_discovery_dir+"/denovo_paths"
                                                for map_with_discovery_dir in input.map_with_discovery_dirs]
    singularity: config["make_prg_dependencies_img"]
    log:
        "logs/add_denovo_paths/{technology}/{coverage}x/{sub_strategy}/{clustering_tool}/{gene}.log"
    script:
        "../scripts/add_denovo_paths.py"



def aggregate_prgs_input(wildcards):
    input_files = []
    for tool, gene in TOOL_MSA_PAIR:
        input_files.append(
            f"analysis/{wildcards.technology}/{wildcards.coverage}x/{wildcards.sub_strategy}/prgs/{tool}/{gene}.prg.fa"
        )

    return input_files


rule aggregate_prgs:
    input:
        prgs = aggregate_prgs_input,
    output:
        "analysis/{technology}/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 500 * attempt
    params:
        original_prg = "data/prgs/ecoli_pangenome_PRG_210619.fa"
    run:
        import fileinput
        with open(output[0], "w") as fout, fileinput.input(input.prgs) as fin:
            for line in fin:
                is_header = line.startswith(">")
                # strip and add newline in case some lines are missing newline
                stripped_line = line.rstrip()
                line_ends_digit = stripped_line[-1].isdigit()
                if (not is_header) and line_ends_digit:
                    stripped_line += " "
                fout.write(stripped_line + "\n")

        # check original prg and new prg have the same number of sequences
        prgs_in_original = 0
        with open(params.original_prg) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_original += 1

        prgs_in_new = 0
        with open(output[0]) as fh:
            for line in fh:
                if line.startswith(">"):
                    prgs_in_new += 1

        assert prgs_in_original == prgs_in_new, f"Original PRG ({prgs_in_original}) and new PRG ({prgs_in_new}) dont have the same number of entries!"
