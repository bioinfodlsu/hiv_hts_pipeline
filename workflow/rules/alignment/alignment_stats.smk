rule vg_index_with_prune:
    input:
        visualize = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome_vg.svg",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.pruned.xg",
        gcsa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.pruned.gcsa",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/",
        pruned_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.pruned.vg",
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/vg_index.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg index -p -x {output.xg} {input.vg} > {log} 2>&1 &&
        vg prune -r {input.vg} > {params.pruned_vg} &&
        vg index -p -g {output.gcsa} {params.pruned_vg} > {log} 2>&1
        """


rule vg_align_stats:
    input:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/alns.gam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/mapped.stats.txt"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/vg_align_stats_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg stats --alignments {input.gam} > {output} 2> {log}
        """

rule gam_to_json:
    input:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/alns.gam",
    output:
        json = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/mapped.json"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
    conda:
        "../../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/gam_to_json_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg view -a {input.gam} -j > {output.json} 2> {log}
        """

rule vg_pack:
    input:
        # json = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/mapped.json",
        gam_stats = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/mapped.stats.txt",
        # xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/alns.gam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.pack",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}",
        node_coverage_out = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.node_coverage.txt",
    conda:
        "../../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/vg_pack_{{sample_id}}.log"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg pack -x {input.vg} -g {input.gam} -o {output} 2> {log}
        """

rule vg_node_coverage:
    input:
        pack = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.pack",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        node_coverage = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.node_coverage.txt"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}",
    conda:
        "../../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/vg_node_coverage_{{sample_id}}.log"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg depth -k {input.pack} {input.vg} > {output.node_coverage} 2> {log}
        """



