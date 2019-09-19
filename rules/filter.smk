rule filtlong:
    input:
        reads = "data/{sample}/{sample}.nanopore.fastq.gz",
        ref = "data/{sample}/{sample}.ref.fa",
    output:
        "data/{sample}/{sample}.{covg}x.nanopore.fastq"
    params:
        mean_q_weight = config['filtlong']['mean_q_weight'],
        min_length = config['filtlong']['min_length'],
    threads: 1
    resources:
        mem_mb = 1000
    singularity: config["filtlong"]["container"]
    log:
        "logs/filtlong/{sample}.{covg}x.log"
    shell:
        """
        bash scripts/downsample_nanopore_reads.sh \
            {input.reads} \
            {input.ref} \
            {wildcards.covg} \
            {output[0]} \
            {params.min_length} \
            {params.mean_q_weight} 2> {log}  
        """
