rule gam_to_hxb2:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns.gam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.gam"
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
        path = "K03455.1",  # HXB2
    conda:
        "../../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/gam_to_hxb2.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/gam_to_hxb2.txt"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&
        vg surject -x {input.vg} -p {params.path} {input.gam} > {output} 2> {log}
        """

rule extract_aln_to_hxb2:
    input:
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns.gam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sam",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}",
        path = "K03455.1",  # HXB2
    conda:
        "../../envs/pggb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/extract_aln_to_hxb2.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/extract_aln_to_hxb2.txt"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&
        vg surject -s -x {input.vg} -p {params.path} {input.gam} > {output} 2> {log}
        """

rule bam_from_surject_sam:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sorted.bam",
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/bam_from_surject_sam.log"
    conda:
        "../../envs/mafft.yaml"
    shell:
        """
        samtools view -S -bh {input} | samtools sort -o {output} 2> {log}
        """

rule index_bam_from_surject:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sorted.bam",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sorted.bai",
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/index_bam_from_surject.log"
    conda:
        "../../envs/mafft.yaml"
    shell:
        """
        samtools index {input} {output} 2> {log}
        """