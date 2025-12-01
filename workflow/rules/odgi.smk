rule odgi_build:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.gfa",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/odgi/pangenome.og",
    params:
        out_dir_odgi = f"{OUT_DIR}/{BATCH_NAME}/odgi",
    conda:
        "../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/odgi/odgi_build.log"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {params.out_dir_odgi} &&
        odgi build -g {input} -o {output} 2> {log}
        """

rule odgi_untangle:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/odgi/pangenome.og",
    output:
        touch(f"{OUT_DIR}/{BATCH_NAME}/odgi/pangenome.untangled.paf"),
    params:
        out_dir_odgi = f"{OUT_DIR}/{BATCH_NAME}/odgi",
    conda:
        "../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/odgi/odgi_untangle.log"
    shell:
        """
        odgi validate -i {input} -P 2> {log} &&\
        odgi untangle -i {input} -p 2> {log}
        """