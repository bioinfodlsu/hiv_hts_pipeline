rule download_refseqs:
    # input: 
    output: 
        touch(f"{OUT_DIR}/{BATCH_NAME}/mafft/download_refseqs.done")
    params:
        email = EMAIL
    log:
        f"logs/{BATCH_NAME}/download_refseqs/download_refseqs.log"
    threads:
        workflow.cores
    conda:
        "../envs/download_refseqs.yaml"
    shell: 
        """
        python workflow/scripts/download_refseqs.py {params.email} > {log} 2>&1
        """

rule clean_refseqs:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/download_refseqs.done"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/cleaned_refs.done"
    conda:
        "../envs/download_refseqs.yaml"
    log:
        f"logs/{BATCH_NAME}/mafft/clean_refseqs.log"
    threads:
        workflow.cores
    shell:
        """
        python workflow/scripts/clean_refseqs.py > {output} 2> {log}
        """

rule fast_tree:
    input:
        lambda wildcards:
            f"{OUT_DIR}/{BATCH_NAME}/mafft/msa_clustered.fasta"
            if config['use_clustered_refs']
            else f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/fast_tree/fast_tree.tre"
    params:
        out_dir = f"{OUT_DIR}/{BATCH_NAME}/fast_tree",
    conda:
        "../envs/fast_tree.yaml"
    log:
        f"logs/{BATCH_NAME}/fast_tree/fast_tree.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir} &&
        FastTree -nt -gtr -gamma {input} > {output} 2> {log}
        """