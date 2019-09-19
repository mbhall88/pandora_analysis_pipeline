from pathlib import Path


checkpoint map_with_discovery:
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
        read_file=$(realpath {input.reads})
        prg_file=$(realpath {input.prg})
        log_file=$(realpath {log})
        mkdir -p {params.outdir}
        cd {params.outdir} || exit 1

        {params.pandora} map --prg_file $prg_file \
            --read_file $read_file \
            --outdir $(pwd) \
            -t {threads} \
            --output_kg \
            --output_covgs \
            --output_vcf \
            --log_level {params.log_level} \
            --genotype \
            --discover > $log_file 2>&1
        """
