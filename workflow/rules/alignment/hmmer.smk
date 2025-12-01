rule refs_nt_to_prot:
    input:
        hxb2 = config['references']['HXB2'],
        refs = f"{OUT_DIR}/{BATCH_NAME}/mafft/concat_refs.fasta",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/concat_refs_prot.faa"
    params:
        out_dir_hmmer = f"{OUT_DIR}/{BATCH_NAME}/hmmer/refs_nt_to_prot",
    conda:
        "../../envs/download_refseqs.yaml"
    log:
        f"logs/{BATCH_NAME}/hmmer/refs_nt_to_prot.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_hmmer} &&
        python workflow/scripts/refs_nt_to_prot.py --hxb2 {input.hxb2} --refs {input.refs} --outdir {params.out_dir_hmmer} 2> {log} &&
        cat {params.out_dir_hmmer}/*_pol_prot.faa > {output} 2>> {log}
        """

rule mafft_prot:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/concat_refs_prot.faa"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/msa_prot.faa"
    params:
        out_dir_hmmer = f"{OUT_DIR}/{BATCH_NAME}/hmmer",
        mafft_opts = "--auto"
    conda:
        "../../envs/mafft.yaml"
    log:
        f"logs/{BATCH_NAME}/hmmer/mafft_prot.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_hmmer} &&
        mafft {params.mafft_opts} {input} > {output} 2> {log}
        """

rule hmmbuild:
    input:
        msa = f"{OUT_DIR}/{BATCH_NAME}/hmmer/msa_prot.faa"
        # msa = f"{OUT_DIR}/{BATCH_NAME}/mafft/msa_trimmed.fasta"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/pol.hmm"
    params:
        out_dir_hmmer = f"{OUT_DIR}/{BATCH_NAME}/hmmer"
    conda:
        "../../envs/hmmer.yaml"
    log:
        f"logs/{BATCH_NAME}/hmmer/hmmbuild.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_hmmer} &&
        hmmbuild {output} {input.msa} 2> {log}
        """

rule make_fastq_unique:
    input:
        reads = lambda wildcards: SAMPLES[wildcards.sample_id]
    output:
        # f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.fastq"
        touch(f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.done")
    params:
        out_dir_fastq = f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}",
    conda:
        "../../envs/download_refseqs.yaml"
    log:
        f"logs/{BATCH_NAME}/fastq/make_fastq_unique_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_fastq} &&
        python workflow/scripts/make_fastq_unique.py --input {input.reads} --outdir {params.out_dir_fastq} --sample_id {wildcards.sample_id} > {log} 2>&1
        """

rule translate_fastq:
    input:
        # reads = lambda wildcards: SAMPLES[wildcards.sample_id],
        done = f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.done",
        reads = f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.fastq"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/{{sample_id}}/reads/{{sample_id}}.6f.faa"
    params:
        out_dir_hmmer = f"{OUT_DIR}/{BATCH_NAME}/hmmer/{{sample_id}}/reads"
    conda:
        "../../envs/hmmer.yaml"
    log:
        f"logs/{BATCH_NAME}/hmmer/translate_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_hmmer} &&
        seqkit translate -f 6 -o {output} {input.reads} 2> {log}
        """

rule hmmsearch:
    input:
        hmm = f"{OUT_DIR}/{BATCH_NAME}/hmmer/pol.hmm",
        read_prot = f"{OUT_DIR}/{BATCH_NAME}/hmmer/{{sample_id}}/reads/{{sample_id}}.6f.faa"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/hmmer/{{sample_id}}/alignments/alns_to_hivmmer.out"
    params:
        out_dir_hmmer = f"{OUT_DIR}/{BATCH_NAME}/hmmer/{{sample_id}}/alignments"
    conda:
        "../../envs/hmmer.yaml"
    log:
        f"logs/{BATCH_NAME}/hmmer/hmmsearch_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_hmmer} &&
        hmmsearch --tblout {params.out_dir_hmmer}/hmmsearch_pol.tbl --cpu {threads} {input.hmm} {input.read_prot} > {output} 2> {log}
        """