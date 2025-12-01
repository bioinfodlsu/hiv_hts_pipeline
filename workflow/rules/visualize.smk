rule visualize_pggb_gfa:
    input:
        rename = f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph/rename.done",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/pangenome_gfa.svg",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg",
        gfa = f"{OUT_DIR}/{BATCH_NAME}/pggb/construct_graph/{BATCH_NAME}_pggb.final.gfa",
    conda:
        "../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/visualize_pggb_gfa.log"
    threads:
        workflow.cores
    shell:
        """
        vg view -dpn -F {params.gfa} | dot -Tsvg -o {output} 2> {log}
        """

rule visualize_gfa_bandage:
    input:
        gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome_bandage.png",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/bandage.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/visualize_gfa_bandage.log"
    threads:
        workflow.cores
    shell:
        """
        Bandage image {input.gfa} {output} --width 2000 --height 1000 2> {log}
        """

rule visualize_vg:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome_vg.svg",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/visualize_vg.log"
    threads:
        workflow.cores
    shell:
        """
        vg view -dp -V {input.vg} | dot -Tsvg -o {output} 2> {log}
        """

rule visualize_vg_gfa:
    input:
        gfa = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome_vg_gfa.svg",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/construct/visualize_vg_gfa.log"
    threads:
        workflow.cores
    shell:
        """
        vg view -dpn -F {input.gfa} | dot -Tsvg -o {output} 2> {log}
        """

rule visualize_gam:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/{{sample_id}}_mapped.gam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/{{sample_id}}_mapped.svg",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map"
    conda:
        "../envs/vg.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/vg_map/visualize_gam_{{sample_id}}.log"
    threads:
        workflow.cores
    shell:
#         vg view -a z.gam | head -1 | vg view -JaG - >first_aln.gam
# vg find -x z.xg -G first_aln.gam | vg view -dA first_aln.gam - | dot -Tpdf -o first_aln.pdf
        """
        mkdir -p {params.out_dir_vg} &&
        vg view -a {input.gam} | head -1 | vg view -JaG - > first_aln.gam &&
        vg find -x {input.vg} -G first_aln.gam | vg view -dA first_aln.gam - | dot -Tsvg -o {output} 2> {log} &&
        rm first_aln.gam
        """