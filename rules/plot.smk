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
