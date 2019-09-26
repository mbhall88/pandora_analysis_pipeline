rule subsample:
    input:
        reads = "data/{sample}/{sample}.nanopore.fastq.gz",
        ref = "data/{sample}/{sample}.ref.fa",
    output:
        "data/{sample}/{sample}.{coverage}x.{sub_strategy}.nanopore.fastq"
    params:
        mean_q_weight = config['subsample']['mean_q_weight'],
        min_length = config['subsample']['min_length'],
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    singularity: config["subsample"]["container"]
    log:
        "logs/subsample/{sub_strategy}/{sample}.{coverage}x.log"
    shell:
        """
        bash scripts/downsample_nanopore_reads.sh \
            {input.reads} \
            {input.ref} \
            {wildcards.coverage} \
            {output[0]} \
            {wildcards.sub_strategy} \
            {params.min_length} \
            {params.mean_q_weight} 2> {log}  
        """
