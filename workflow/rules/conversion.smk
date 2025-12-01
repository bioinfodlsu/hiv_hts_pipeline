rule vg_to_gfa:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/vg_to_gfa.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_vg} &&
        vg view {input.vg} > {output} 2> {log}
        """

# rule vg_to_xg:
#     input:
#         vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
#     output:
#         xg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.xg",
#     params:
#         out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
#     conda:
#         "../envs/pggb.yaml"
#     log:
#         f"logs/{BATCH_NAME}/vg/construct/vg_to_xg.log"
#     threads:
#         workflow.cores
#     shell:
#         """
#         mkdir -p {params.out_dir_vg} &&
#         vg convert -x {input.vg} > {output.xg} 2> {log}
#         """