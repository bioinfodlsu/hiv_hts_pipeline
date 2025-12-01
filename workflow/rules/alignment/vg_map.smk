# rule vg_index:
#     input:
#         visualize = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome_vg.svg",
#         stats = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.stats.txt",
#         vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
#         gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
#     output:
#         xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
#         gcsa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gcsa",
#     params:
#         out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/",
#     conda:
#         "../../envs/vg.yaml"
#     log:
#         f"logs/{BATCH_NAME}/vg/construct/vg_index.log"
#     threads:
#         workflow.cores
#     shell:
#         """
#         mkdir -p {params.out_dir_vg} &&
#         vg index -p -x {output.xg} -g {output.gcsa} {input.vg} >> {log} 2>&1
#         """
rule vg_autoindex_map:
    input:
        gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    output:
        touch(f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/vg_autoindex.done"),
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map",
        workflow = "map",
        prefix = "pangenome-map"
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map/vg_autoindex.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg autoindex --workflow {params.workflow} -g {input.gfa} -p {params.out_dir_vg}/{params.prefix} > {log} 2>&1
        """

# ruleorder: vg_index > vg_index_with_prune

rule vg_map:
    input:
        # done = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/vg_index.done",
        # xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
        # gcsa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gcsa",
        autoindex_done = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/vg_autoindex.done",
        # reads = lambda wildcards: SAMPLES[wildcards.sample_id],
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
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/{{sample_id}}/alns.gam",
    params:
        index_dir = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map",
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/{{sample_id}}",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/vg_map/{{sample_id}}",
        prefix = "pangenome-map"
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map/{{sample_id}}/vg_map.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/vg_map/{{sample_id}}/vg_map.txt"
    threads:
        workflow.cores
    shell:
        r"""
        set -euo pipefail;

        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&

        if [ $(echo "{input.query}" | wc -w) -eq 2 ]; then
            vg map -f {input.read1} -f {input.read2} -x {params.index_dir}/{params.prefix}.xg -g {params.index_dir}/{params.prefix}.gcsa > {output.gam} 2> {log} 
        else
            vg map -f {input.read1} -x {params.index_dir}/{params.prefix}.xg -g {params.index_dir}/{params.prefix}.gcsa > {output.gam} 2> {log} 
        fi
        """