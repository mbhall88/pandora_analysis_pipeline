rule minimap2:
    input:
        target = "data/{sample}/{sample}.ref.fa",
        query = "data/{sample}/{sample}.{covg}x.nanopore.fastq"
    output:
        "analysis/{covg}x/alignments/{sample}.sorted.bam"
    log:
        "logs/minimap2/{sample}.{covg}x.log"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000
    shell:
        """
        minimap2 -t {threads} \
            -ax map-ont \
            {input.target} {input.query} | \
                samtools sort -@ {threads} -o {output} - 2> {log}
        """

rule minimap2_original_nanopore_data:
    input:
        target = "data/{sample}/{sample}.ref.fa",
        query = "data/{sample}/{sample}.nanopore.fastq.gz"
    output:
        "analysis/alignments/{sample}.sorted.bam"
    log:
        "logs/minimap2/{sample}.log"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000
    shell:
        """
        minimap2 -t {threads} \
            -ax map-ont \
            {input.target} {input.query} | \
                samtools sort -@ {threads} -o {output} - 2> {log}
        """

rule plot_filtered_data:
    input:
        "analysis/{covg}x/alignments/{sample}.sorted.bam"
    output:
        "analysis/{covg}x/plots/{sample}/NanoPlot-report.html"
    singularity: config["plot"]["container"]
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000
    params:
        (
            " --loglength"
            " --verbose"
        )
    log:
        "logs/plot/{sample}.{covg}x.log"
    shell:
        "NanoPlot --threads {threads} {params} --bam {input[0]} --outdir $(dirname {output[0]}) 2> {log}"

rule plot_original_data:
    input:
        "analysis/alignments/{sample}.sorted.bam"
    output:
        "analysis/plots/{sample}/NanoPlot-report.html"
    singularity: config["plot"]["container"]
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000
    params:
        (
            " --loglength"
            " --verbose"
        )
    log:
        "logs/plot/{sample}.log"
    shell:
        "NanoPlot --threads {threads} {params} --bam {input[0]} --outdir $(dirname {output[0]}) 2> {log}"


