rule index_original_prg:
    input:
        "data/prgs/ecoli_pangenome_PRG_210619.fa"
    output:
        index = "data/prgs/ecoli_pangenome_PRG_210619.fa.k15.w14.idx",
        kmer_prgs = directory("data/prgs/kmer_prgs")
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000
    params:
        pandora = "/nfs/research1/zi/mbhall/Software/pandora/build-release/pandora"
    log:
        "logs/index_original_prg.log"
    shell:
        "{params.pandora} index -t {threads} {input} > {log} 2>&1"

rule index_prg_updated_with_denovo_paths:
    input:
        "analysis/{coverage}x/prgs/denovo_updated.prg.fa"
    output:
        index = "analysis/{coverage}x/prgs/denovo_updated.prg.fa.k15.w14.idx",
        kmer_prgs = directory("analysis/{coverage}x/prgs/kmer_prgs")
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000
    params:
        pandora = "/nfs/research1/zi/mbhall/Software/pandora/build-release/pandora"
    log:
        "logs/{coverage}x/prgs/index_prg_updated_with_denovo_paths.log" 
    shell:
        "{params.pandora} index -t {threads} {input} > {log} 2>&1"
