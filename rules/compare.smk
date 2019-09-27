from pathlib import Path

rule create_tsv_for_reads:
    input:
        expand("data/{sample}/{sample}.{{coverage}}x.{{sub_strategy}}.nanopore.fastq",
               sample=config["samples"])
    output:
        tsv="data/samples.{coverage}x.{sub_strategy}.tsv"
    threads: 1
    resources:
        mem_mb=200
    log:
        "logs/create_tsv_for_reads/{coverage}x/{sub_strategy}.log"
    shell:
        """
        for path in {input}
        do
            filename=$(basename $path)
            sample_name=${{filename/.fastq}}
            echo -e \"$sample_name\t$(realpath $path)\" >> {output.tsv} 2>> {log}
        done
        """

rule compare_with_denovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg="analysis/{coverage}x/{sub_strategy}/prgs/denovo_updated.prg.fa",
        prg_index=rules.index_prg_updated_with_denovo_paths.output.index,
    output:
        vcf="analysis/{coverage}x/{sub_strategy}/compare_with_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref="analysis/{coverage}x/{sub_strategy}/compare_with_denovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        pandora=config["pandora_executable"],
        log_level="debug",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
    log:
        "logs/compare_with_denovo/{coverage}x/{sub_strategy}.log"
    shell:
        """
        {params.pandora} compare --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype \
            --log_level {params.log_level} > {log} 2>&1
        """

rule compare_no_denovo:
    input:
        read_index=rules.create_tsv_for_reads.output.tsv,
        prg=config["original_prg"],
        prg_index=rules.index_original_prg.output.index,
    output:
        vcf="analysis/{coverage}x/{sub_strategy}/compare_no_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref="analysis/{coverage}x/{sub_strategy}/compare_no_denovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000
    params:
        pandora=config["pandora_executable"],
        log_level="debug",
        outdir=lambda wildcards, output: str(Path(output.vcf).parent),
    log:
        "logs/compare_no_denovo/{coverage}x/{sub_strategy}.log"
    shell:
        """
        {params.pandora} compare \
            --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir {params.outdir} \
            -t {threads} \
            --genotype \
            --log_level {params.log_level} > {log} 2>&1
        """
