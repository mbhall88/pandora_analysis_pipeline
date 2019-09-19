rule create_tsv_for_reads:
    input:
        expand("data/{sample}/{sample}.{{coverage}}x.nanopore.fastq", sample=config["samples"])
    output:
        "data/samples.{coverage}x.tsv"
    threads: 1
    resources:
        mem_mb = 200
    log:
        "logs/create_tsv_for_reads.{coverage}x.log"
    shell:
        """
        for path in {input}
        do
            filename=$(basename $path)
            sample_name=${{filename/.fastq}}
            echo -e \"$sample_name\t$(realpath $path)\" >> {output} 2>> {log}
        done
        """ 

rule compare_with_denovo:
    input:
        read_index = "data/samples.{coverage}x.tsv",
        prg = "analysis/{coverage}x/prgs/denovo_updated.prg.fa",
        prg_index = "analysis/{coverage}x/prgs/denovo_updated.prg.fa.k15.w14.idx",
    output:
        vcf = "analysis/{coverage}x/compare_with_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref = "analysis/{coverage}x/compare_with_denovo/pandora_multisample.vcf_ref.fa"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000
    params:
        pandora = "/nfs/research1/zi/mbhall/Software/pandora/build-release/pandora"
    log:
        "data/{coverage}x/compare_with_denovo.log"
    shell:
        """
        outdir=$(dirname {output.vcf})
        {params.pandora} compare --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir $outdir \
            -t {threads} \
            --genotype \
            --log_level debug > {log} 2>&1
        """

rule compare_without_denovo:
    input:
        read_index = "data/samples.{coverage}x.tsv",
        prg = "data/prgs/ecoli_pangenome_PRG_210619.fa",
        prg_index = "data/prgs/ecoli_pangenome_PRG_210619.fa.k15.w14.idx",
    output:
        vcf = "analysis/{coverage}x/compare_without_denovo/pandora_multisample_genotyped.vcf",
        vcf_ref = "analysis/{coverage}x/compare_without_denovo/pandora_multisample.vcf_ref.fa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000
    params:
        pandora = "/nfs/research1/zi/mbhall/Software/pandora/build-release/pandora"
    log:
        "data/{coverage}x/compare_without_denovo.log"
    shell:
        """
        outdir=$(dirname {output.vcf})
        {params.pandora} compare --prg_file {input.prg} \
            --read_index {input.read_index} \
            --outdir $outdir \
            -t {threads} \
            --genotype \
            --log_level debug > {log} 2>&1
        """

