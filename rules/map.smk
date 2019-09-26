from pathlib import Path


rule map_with_discovery:
    input:
        prg = config["original_prg"],
        index = rules.index_original_prg.output.index,
        reads = "data/{sample}/{sample}.{coverage}x.{sub_strategy}.nanopore.fastq",
    output:
        denovo_dir = directory("analysis/{coverage}x/{sub_strategy}/{sample}/map_with_discovery/denovo_paths"),
        consensus = "analysis/{coverage}x/{sub_strategy}/{sample}/map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{coverage}x/{sub_strategy}/{sample}/map_with_discovery/pandora_genotyped.vcf",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000
    params:
        outdir = lambda wildcards, output: str(Path(output.consensus).parent),
        pandora = config["pandora_executable"],
        log_level = "debug",
    log:
        "logs/map_with_discovery/{coverage}x/{sub_strategy}/{sample}.log"
    shell:
        """
        {params.pandora} map --prg_file {input.prg} \
            --read_file {input.reads} \
            --outdir {params.outdir} \
            -t {threads} \
            --output_kg \
            --output_covgs \
            --output_vcf \
            --log_level {params.log_level} \
            --genotype \
            --discover > {log} 2>&1
        """
