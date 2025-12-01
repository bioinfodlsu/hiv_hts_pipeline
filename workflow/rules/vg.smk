rule vg_construct_msa:
    input:
        # msa = f"{OUT_DIR}/{BATCH_NAME}/mafft/msa_trimmed.fasta"
        msa = f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/vg_construct_msa.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/construct/vg_construct_msa.txt"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&
        vg construct -M {input.msa} -p > {output} 2> {log}
        """

rule vg_prune:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.pruned.vg",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/vg_prune.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg prune {input.vg} > {output} 2> {log}
        """

rule vg_stats:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.stats.txt",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct",
        stats = "--size --length --self-loops --subgraphs --heads --tails --nondeterm --is-acyclic"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/vg_stats.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg stats {params.stats} {input.vg} > {output} 2> {log}
        """