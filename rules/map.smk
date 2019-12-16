rule map_with_discovery:
    input:
        prg=config["original_prg"],
        index=rules.index_original_prg.output.index,
        reads="data/{sample}/{sample}.{coverage}x.{sub_strategy}.{technology}.fastq",
        ref="data/{sample}/{sample}.ref.fa",
    output:
        outdir=directory("analysis/{technology}/{coverage}x/{sub_strategy}/{sample}/map_with_discovery")
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 30000,
    params:
        pandora=config["pandora_executable"],
        log_level="debug",
        use_discover=True,
    log:
        "logs/map_with_discovery/{technology}/{coverage}x/{sub_strategy}/{sample}.log"
    shell:
        """
        bash scripts/pandora_map.sh {params.pandora} {input.prg} \
            {input.reads} \
            {output.outdir} \
            {threads} \
            {wildcards.technology} \
            {params.log_level} \
            {log} \
            {params.use_discover} \
            {input.ref}
        """
