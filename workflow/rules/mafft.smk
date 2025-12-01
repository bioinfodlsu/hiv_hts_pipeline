rule cat_fasta:
    input:
        refs = lambda wildcards: REFERENCES if config['use_own_references'] 
        # else f"{OUT_DIR}/{BATCH_NAME}/mafft/download_refseqs.done"
        else f"{OUT_DIR}/{BATCH_NAME}/mafft/cleaned_refs.done"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta"
    params:
        in_dir = f"test_data/refseqs_cleaned",
        out_dir_mafft = f"{OUT_DIR}/{BATCH_NAME}/mafft"
    conda:
        "../envs/mafft.yaml"
    log:
        f"logs/{BATCH_NAME}/mafft/cat_fasta.log"
    threads:
        workflow.cores
    shell:
        (
            f"mkdir -p {{params.out_dir_mafft}} && "
            f"cat {{input.refs}} > {{output}} 2> {{log}}"
            if config['use_own_references']
            else
            f"mkdir -p {{params.out_dir_mafft}} && "
            f"cat {{params.in_dir}}/K03455.fasta "
            f"$(ls {{params.in_dir}}/*.fasta | grep -v 'K03455.fasta') "
            f"> {{output}} 2> {{log}}"
        )

rule mafft_msa:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
    params:
        out_dir_mafft = f"{OUT_DIR}/{BATCH_NAME}/mafft",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/mafft",
        mafft_opts = "--auto"
    conda:
        "../envs/mafft.yaml"
    log:
        f"logs/{BATCH_NAME}/mafft/mafft_msa.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/mafft/mafft_msa.txt"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_mafft} {params.benchmark_dir} &&
        mafft {input} > {output} 2> {log}
        """

rule msa_trim_to_pol:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/msa_trimmed.fasta"
    params:
        out_dir_mafft = f"{OUT_DIR}/{BATCH_NAME}/mafft",
        offset = 0
    conda:
        "../envs/download_refseqs.yaml"
    log:
        f"logs/{BATCH_NAME}/mafft/msa_trim.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_mafft} &&
        python workflow/scripts/msa_trim_to_pol.py --input {input} --offset {params.offset} --output {output} 2> {log}
        """

rule msa_index:
    input:
        # f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fasta"
        lambda wildcards:
            f"{OUT_DIR}/{BATCH_NAME}/mafft/clustered_refs.fasta"
            if config['use_clustered_refs']
            # else f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta"
            else f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/mafft/msa.fai"
    params:
        out_dir_mafft = f"{OUT_DIR}/{BATCH_NAME}/mafft"
    conda:
        "../envs/mafft.yaml"
    log:
        f"logs/{BATCH_NAME}/mafft/msa_index.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_mafft} &&
        samtools faidx {input} > {output} 2> {log}
        """