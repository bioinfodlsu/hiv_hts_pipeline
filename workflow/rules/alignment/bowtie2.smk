rule index_bowtie2:
    input:
        reference=f"test_data/refseqs_cleaned/{{reference_name}}.fasta",
    output:
        touch(
            f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/index/{{reference_name}}_index.done"
        ),
    params:
        index_basename=f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/index/{{reference_name}}_index",
    log:
        f"logs/{BATCH_NAME}/bowtie2/{{sample_id}}/index/{{reference_name}}_index.log",
    conda:
        "../../envs/linear.yaml"
    shell:
        "bowtie2-build {input.reference} {params.index_basename} > {log} 2>&1"


rule align_bowtie2:
    input:
        reference_flag=f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/index/{{reference_name}}_index.done",
        query = lambda wildcards: config["samples"][wildcards.sample_id],
        read1=(
            f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.done"
            if PREPROCESS_READS
            else lambda wildcards: config["samples"][wildcards.sample_id][0]
        ),
        read2=(
            lambda wildcards: config["samples"][wildcards.sample_id][1]
            if is_paired(wildcards.sample_id)
            else lambda wildcards: config["samples"][wildcards.sample_id][0]
        )
    output:
        f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/alignments/alns_to_{{reference_name}}.sam",
    params:
        out_dir_bowtie2 = f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/alignments",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/bowtie2/{{sample_id}}/alignments",
        index_basename=rules.index_bowtie2.params.index_basename,
        preset="sensitive-local",
        read1=(f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.fastq"
                if PREPROCESS_READS
                else lambda wildcards: config["samples"][wildcards.sample_id][0]),
    log:
        f"logs/{BATCH_NAME}/bowtie2/{{sample_id}}/alignments/alns_to_{{reference_name}}.log",
    benchmark:
        f"benchmarks/{BATCH_NAME}/bowtie2/{{sample_id}}/alignments/alns_to_{{reference_name}}.txt",
    conda:
        "../../envs/linear.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_dir_bowtie2} {params.benchmark_dir} &&

        if [ $(echo "{input.query}" | wc -w) -eq 2 ]; then
            # Paired-end
            bowtie2 --local --{params.preset} -x {params.index_basename} -1 {input.read1} -2 {input.read2} -S {output} > {log} 2>&1
        else
            # Single-end
            bowtie2 --local --{params.preset} -x {params.index_basename} -U {input.read1} -S {output} > {log} 2>&1
        fi
        """

rule convertSAM:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/alignments/alns_to_{{reference_name}}.sam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/bams/alns_to_{{reference_name}}.bam",
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/bams/alns_to_{{reference_name}}.log",
    threads: workflow.cores
    conda:
        "../../envs/mafft.yaml"
    shell:
        """
        samtools view -S -bh {input} > {output} 2> {log}
        """


rule sortBAM:
    input:
        # bai = f"{OUT_DIR}/{BATCH_NAME}/bowtie2/{{sample_id}}/bam/alns_to_{{reference_name}}.bai",
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/bams/alns_to_{{reference_name}}.bam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bam",
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bam.log",
    conda:
        "../../envs/mafft.yaml"
    shell:
        "samtools sort {input} -o {output} 2> {log}"


rule indexBAM:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bai",
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bai.log",
    conda:
        "../../envs/mafft.yaml"
    shell:
        "samtools index {input} > {output} 2> {log}"
