rule vg_sim:
    input:
        xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/alignment_simulated.gam"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim",
        num_reads = 100000,  # Number of reads to simulate
        read_length = 100,  # Length of each read
        substitution_rate = 0.1,  # Substitution rate
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_sim/vg_sim.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg sim -x {input.xg} -n {params.num_reads} -l {params.read_length} -e {params.substitution_rate} -a > {output} 2> {log} 
        """

rule extract_simulated_reads:
    input:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/alignment_simulated.gam"
    output:
        reads = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/simulated_reads.fastq.gz"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim",
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_sim/extract_simulated_reads.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg view -X {input.gam} | gzip > {output} 2> {log}
        """

rule vg_map_compare:
    input:
        sim_gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/alignment_simulated.gam",
        xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
        gcsa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gcsa",
        reads = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/simulated_reads.fastq.gz",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}.gam"
    output:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/{{sample_id}}_mapped.gam",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim",
    log:
        f"logs/{BATCH_NAME}/vg/vg_sim/vg_map_compare_{{sample_id}}.log"
    conda:
        "../envs/vg.yaml"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg map --compare {input.sim_gam} -G {input.gam} -x {input.xg} -g {input.gcsa} > {output.gam} 2> {log}
        """

rule vg_map_eval:
    input:
        # done = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/vg_index.done",
        xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
        gcsa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gcsa",
        reads = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/simulated_reads.fastq.gz",
    output:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}.gam"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval",
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map_eval/vg_map_eval_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg map -f {input.reads} -x {input.xg} -g {input.gcsa} > {output.gam} 2> {log} 
        """

rule vg_eval_align_stats:
    input:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}.gam"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}_mapped.stats.txt"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval",
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map_eval/vg_align_stats_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg stats --alignments {input.gam} > {output} 2> {log}
        """

rule eval_to_json:
    input:
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}.gam"
    output:
        json = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}_mapped.json"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval",
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map_eval/gam_to_json_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg view -a {input.gam} -j > {output.json} 2> {log}
        """

rule vg_pack_eval:
    input:
        json = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}_mapped.json",
        gam_stats = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}_mapped.stats.txt",
        xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map_eval/{{sample_id}}.gam"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval/{{sample_id}}_mapped.pack"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval",
        node_coverage_out = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval/{{sample_id}}_mapped.node_coverage.txt",
    conda:
        "../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_pack_eval/vg_pack_{{sample_id}}.log"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg pack -x {input.vg} -g {input.gam} -o {output} 2> {log}
        """

rule vg_node_coverage_eval:
    input:
        pack = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval/{{sample_id}}_mapped.pack",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        compare_gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_sim/{{sample_id}}_mapped.gam",
    output:
        node_coverage = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval/{{sample_id}}_mapped.node_coverage.txt"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_pack_eval",
    conda:
        "../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_pack_eval/vg_node_coverage_{{sample_id}}.log"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg depth -k {input.pack} {input.vg} > {output.node_coverage} 2> {log}
        """