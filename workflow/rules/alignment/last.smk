from os import path
import itertools


def getMean(sample_id, reference_name, out_dir, batch_name):
    #     est_path = "{}/last_frag_stat/{}.frag_len_est".format(config["out_dir"],wildcards.sample_id)
    est_path = f"{out_dir}/{batch_name}/last/{sample_id}/last_frag_stat/{sample_id}_to_{reference_name}.frag_len_est"
    if path.exists(est_path):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated mean distance"):
                    return line.rstrip().split()[-1]
    else:
        return -1


def getSTD(sample_id, reference_name, out_dir, batch_name):
    est_path = f"{out_dir}/{batch_name}/last/{sample_id}/last_frag_stat/{sample_id}_to_{reference_name}.frag_len_est"
    if path.exists(est_path):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated standard deviation of distance"):
                    return line.rstrip().split()[-1]
    else:
        return -1


rule last_db:
    input:
        reference=f"test_data/refseqs_cleaned/{{reference_name}}.fasta",
    output:
        touch(
            f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.done"
        ),
    params:
        index_basename=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index",
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.log",
    conda:
        "../../envs/linear.yaml"
    shell:
        "lastdb {params.index_basename} {input.reference} > {log} 2>&1"


rule last_frag_stat_est:
    input:
        reference_flag=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.done",
        query=lambda wildcards: (
            f"{OUT_DIR}/{BATCH_NAME}/fastq/{wildcards.sample_id}/{wildcards.sample_id}_unique.done"
            if PREPROCESS_READS
            else config["samples"][wildcards.sample_id]
        ),
    params:
        outdir=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_frag_stat",
        index_basename=rules.last_db.params.index_basename,
        last_al=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_frag_stat/{{sample_id}}_to_{{reference_name}}.frag_len_est",
    output:
        frag_len_est=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_frag_stat/{{sample_id}}_to_{{reference_name}}.frag_len_est",
    conda:
        "../../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/last_frag_stat/{{sample_id}}_to_{{reference_name}}_frag_len_est.log",
    threads: workflow.cores
    shell:
        """
        set +o pipefail;
        mkfifo {params.outdir}/{wildcards.sample_id}_1 |
        mkfifo {params.outdir}/{wildcards.sample_id}_2 |

        gzip -cdf {input.query[0]} | head -n 400000 > {params.outdir}/{wildcards.sample_id}_1 &
        gzip -cdf {input.query[1]} | head -n 400000 > {params.outdir}/{wildcards.sample_id}_2 &

        fastq-interleave {params.outdir}/{wildcards.sample_id}_1 {params.outdir}/{wildcards.sample_id}_2 | lastal -Q1 -i1 {params.index_basename} | last-pair-probs -e > {output.frag_len_est} 2> {log}

        rm {params.outdir}/{wildcards.sample_id}_1 {params.outdir}/{wildcards.sample_id}_2
        """


rule create_dict:
    input:
        reference=f"test_data/refseqs_cleaned/{{reference_name}}.fasta",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.dict",
    conda:
        "../../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index_dict.log",
    shell:
        "picard CreateSequenceDictionary R={input.reference} O={output} > {log} 2>&1"


rule last_score_training:
    input:
        reference_flag=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.done",
        query=(
            f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.done"
            if PREPROCESS_READS
            else lambda wildcards: config["samples"][wildcards.sample_id]
        ),
    params:
        outdir=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_training",
        benchmark_dir=f"benchmarks/{BATCH_NAME}/last/{{sample_id}}/last_training",
        index_basename=rules.last_db.params.index_basename,
        last_sample=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_training/{{sample_id}}_to_{{reference_name}}_sample.fastq",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_training/{{sample_id}}_to_{{reference_name}}.scoring_scheme",
    conda:
        "../../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/last_training/{{sample_id}}_to_{{reference_name}}_scoring_scheme.log",
    benchmark:
        f"benchmarks/{BATCH_NAME}/last/{{sample_id}}/last_training/{{sample_id}}_to_{{reference_name}}_scoring_scheme.txt",
    threads: workflow.cores
    shell:
        """
        set +o pipefail;
        mkdir -p {params.outdir} {params.benchmark_dir} &&
        gzip -cdf {input.query} | head -n 400000 > {params.last_sample}
        last-train --verbose --revsym --matsym --gapsym --sample-number=400000 -Q1  {params.index_basename} {params.last_sample} > {output}
        """

rule last_align:
    input:
        reference_flag=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.done",
        query=(
            f"{OUT_DIR}/{BATCH_NAME}/fastq/{{sample_id}}/{{sample_id}}_unique.done"
            if PREPROCESS_READS
            else lambda wildcards: config["samples"][wildcards.sample_id]
        ),
        scoring_scheme=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/last_training/{{sample_id}}_to_{{reference_name}}.scoring_scheme",
        # Only require frag_len_est for paired-end samples
        frag_len_est=lambda wildcards: (
            f"{OUT_DIR}/{BATCH_NAME}/last/{wildcards.sample_id}/last_frag_stat/{wildcards.sample_id}_to_{wildcards.reference_name}.frag_len_est"
            if is_paired(wildcards.sample_id)
            else []
        ),
    params:
        outdir=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/alignments",
        index_basename=rules.last_db.params.index_basename,
        # Only get mean/std for paired-end samples
        mean=lambda wildcards: (
            getMean(wildcards.sample_id, wildcards.reference_name, OUT_DIR, BATCH_NAME)
            if is_paired(wildcards.sample_id)
            else None
        ),
        std=lambda wildcards: (
            getSTD(wildcards.sample_id, wildcards.reference_name, OUT_DIR, BATCH_NAME)
            if is_paired(wildcards.sample_id)
            else None
        ),
        benchmark_dir=f"benchmarks/{BATCH_NAME}/last/{{sample_id}}/alignments",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}.maf",
    conda:
        "../../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}.log",
    benchmark:
        f"benchmarks/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}.txt",
    threads: workflow.cores
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir} {params.benchmark_dir} &&

        if [ $(echo "{input.query}" | wc -w) -eq 2 ]; then
            # Paired-end alignments
                fastq-interleave {input.query} |
                lastal -j4 -Q1 -i1 -p {input.scoring_scheme} {params.index_basename} |
                last-pair-probs -f {params.mean} -s {params.std} -d 0.1 > {output} 2> {log}
        else
                lastal -j4 -Q1 -i1 -p {input.scoring_scheme} {params.index_basename} {input.query} |
                last-map-probs > {output} 2> {log}
        fi
        """

rule maf_to_sam:
    input:
        maf=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}.maf",
        ref_dict=f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/index/{{reference_name}}_index.dict",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}.sam",
    conda:
        "../../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/last/{{sample_id}}/alignments/alns_to_{{reference_name}}_maf_to_sam.log",
    shell:
        "maf-convert sam -f {input.ref_dict} {input.maf} > {output} 2> {log}"