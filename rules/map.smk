from pathlib import Path


checkpoint map_with_discovery:
    input:
        prg = "data/prgs/ecoli_pangenome_PRG_210619.fa",
        index = "data/prgs/ecoli_pangenome_PRG_210619.fa.k15.w14.idx",
        reads = "data/{sample}/{sample}.{covg}x.{sub_strategy}.nanopore.fastq",
    output:
        denovo_dir = directory("analysis/{covg}x/{sub_strategy}/{sample}/map_with_discovery/denovo_paths"),
        consensus = "analysis/{covg}x/{sub_strategy}/{sample}/map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{covg}x/{sub_strategy}/{sample}/map_with_discovery/pandora_genotyped.vcf",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000
    params:
        outdir = lambda wildcards, output: str(Path(output.consensus).parent),
        pandora = config["pandora_executable"],
    log:
        "logs/map_with_discovery/{covg}x/{sub_strategy}/{sample}.log"
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
            --log_level debug \
            --genotype \
            --discover > $log_file 2>&1
        """
