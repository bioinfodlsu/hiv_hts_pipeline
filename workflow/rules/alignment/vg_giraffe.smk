rule vg_autoindex:
    input:
        gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    output:
        touch(f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/vg_autoindex.done"),
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe",
        workflow = "giraffe",
        prefix = "pangenome-giraffe"
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_giraffe/vg_autoindex.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg autoindex --workflow {params.workflow} -g {input.gfa} -p {params.out_dir_vg}/{params.prefix} > {log} 2>&1
        """

rule vg_giraffe:
    input:
        autoindex_done = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/vg_autoindex.done",
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
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/alns.gam",
    params:
        index_dir = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe",
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}",
        prefix = "pangenome-giraffe",
        preset = "default",
        rescue = lambda wildcards: "none" if is_paired(wildcards.sample_id) else "dozeu",  # dozeu, none
        # mean = 600.0,
        # sd = 100.0,
        # verbose = "--show-work --track-provenance"
        verbose = ""
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/vg_giraffe.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/vg_giraffe.txt"
    threads:
        workflow.cores
    shell:
        r"""
        set -euo pipefail;

        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&

        if [ $(echo "{input.query}" | wc -w) -eq 2 ]; then
            vg giraffe {params.verbose} --progress --log-reads --threads {threads} -A {params.rescue} --parameter-preset {params.preset} -Z {params.index_dir}/{params.prefix}.giraffe.gbz -m {params.index_dir}/{params.prefix}.shortread.withzip.min -z {params.index_dir}/{params.prefix}.shortread.zipcodes -d {params.index_dir}/{params.prefix}.dist -f {input.read1} -f {input.read2} > {output.gam} 2> {log}
        else
            vg giraffe {params.verbose} --progress --log-reads --threads {threads} -A {params.rescue} --parameter-preset {params.preset} -Z {params.index_dir}/{params.prefix}.giraffe.gbz -m {params.index_dir}/{params.prefix}.shortread.withzip.min -z {params.index_dir}/{params.prefix}.shortread.zipcodes -d {params.index_dir}/{params.prefix}.dist -f {input.read1} > {output.gam} 2> {log}
        fi
        """