rule pggb_construct_graph:
    input:
        index = f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fai",
        msa = lambda wildcards: f"{OUT_DIR}/{BATCH_NAME}/clustering/clustered_refs.fasta" if config['use_clustered_refs'] else f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta",
        mafft = f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
    output:
        touch(f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph/construct_graph.done"),
        # and the main GFA file that pggb will generate inside it
        # gfa = f"{OUT_DIR}/{BATCH_NAME}/pggb/pangenome.gfa",
    params:
        dir = directory(f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph"),
        p = 98,                 # min avg. nucleotide identity
        s = 256,                # segment length
        n = len(REFERENCES),    # number of haplotypes
        k = 500,
        l = 500,

    conda:
        "../envs/pggb.yaml",
    log:
        f"logs/{BATCH_NAME}/pggb/pggb_construct_graph.log",
    threads:
        workflow.cores,
    shell:
        """
        rm -rf {params.dir} &&
        mkdir -p {params.dir} &&
        pggb \
            -i {input.msa} \
            -o {params.dir} \
            -n {params.n} \
            -t {threads} \
            -p {params.p} \
            -s {params.s} \
        > {log} 2>&1
        """
        
rule pggb_rename:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph/construct_graph.done",
    output:
        touch(f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph/rename.done"),
    params:
        dir = directory(f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph"),
        new_name = f"{BATCH_NAME}_pggb"
    log:
        f"logs/{BATCH_NAME}/pggb/pggb_rename.log",
    threads:
        workflow.cores,
    conda:
        "../envs/pggb.yaml",
    shell:
        """
        python workflow/scripts/rename_pggb_outputs.py {params.dir} {params.new_name} > {log} 2>&1
        """
