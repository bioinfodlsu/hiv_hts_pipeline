rule vg_call:
    input:
        # svg = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_map/{{sample_id}}_mapped.svg",
        node_coverage = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.node_coverage.txt",
        pack = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/vg_pack/{{sample_id}}/mapped.pack",
        vg = f"{OUT_DIR}/{BATCH_NAME}/vg/construct/pangenome.vg",
        gam = f"{OUT_DIR}/{BATCH_NAME}/vg/vg_giraffe/{{sample_id}}/alns.gam",
    output:
        vcf = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/vg_call.vcf",
    params:
        out_dir_vg = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants",
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/vg_call.log",
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/vg_call.txt",
    conda:
        "../envs/pggb.yaml"
    threads:
        workflow.cores/len(SAMPLES)
    shell:
        """
        mkdir -p {params.out_dir_vg} {params.benchmark_dir} &&
        vg call {input.vg} -k {input.pack} > {output.vcf} 2> {log}
        """

rule lofreq_linear:
    input:
        reference=f"test_data/refseqs_cleaned/{{reference_name}}.fasta",
        bam = f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bam",
        bai = f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/sorted/alns_to_{{reference_name}}.sorted.bai",
    output:
        vcf=f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.vcf",
    params:
        outdir=f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants",
        benchmark_dir = f"benchmarks/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants",
    conda:
        "../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.log",
    benchmark:
        f"benchmarks/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.txt",
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir} {params.benchmark_dir} &&
        lofreq faidx {input.reference}
        lofreq call-parallel --no-default-filter --pp-threads {threads} -f {input.reference} -o {output} {input.bam} 2> {log}
        """

rule lofreq_vg_sam:
    input:
        reference = "test_data/references/hxb2.fasta",
        bam = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sorted.bam",
        bai = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/alns_to_HXB2.sorted.bai",
    output:
        vcf = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq.vcf",
    params:
        outdir=f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants",
        benchmark_dir = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/",
    conda:
        "../envs/linear.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq_vg_sam.log",
    benchmark:
        f"benchmarks/{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq_vg_sam.txt",
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir} {params.benchmark_dir} &&
        lofreq faidx {input.reference}
        lofreq call-parallel --no-default-filter --pp-threads {threads} -f {input.reference} -o {output} {input.bam} 2> {log}
        """

rule vcf_to_aavf_linear:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.vcf"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.aavf" 
    params:
        out_dir = f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants"
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/vcf_to_aavf_alns_to_{{reference_name}}.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/vcf_to_aavf_alns_to_{{reference_name}}.txt"
    shell:
        """
        mkdir -p {params.out_dir} &&
        python workflow/scripts/vcf_to_aavf_2.py --vcf {input} --aavf {output} 2> {log}
        """

rule vcf_to_aavf_vg:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq.vcf",
    output:
       f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq.aavf",
    params:
        out_dir = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants",
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/vcf_to_aavf_vg.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/vcf_to_aavf_vg.txt"
    shell:
        """
        mkdir -p {params.out_dir} &&
        python workflow/scripts/vcf_to_aavf_2.py --vcf {input} --aavf {output} 2> {log}
        """

rule aavf_filtered_linear:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.aavf"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/alns_to_{{reference_name}}_{{lofreq_value}}.aavf"
    params:
        out_dir=lambda wc: f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{wc.sample_id}/variants_filtered/{wc.lofreq_group}",
        flag=lambda wc: "--min-acf" if wc.lofreq_group == "AF" else "--min-acc",
        value=lambda wc: wc.lofreq_value
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/alns_to_{{reference_name}}_{{lofreq_value}}.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/alns_to_{{reference_name}}_{{lofreq_value}}.txt"
    shell:
        """
        mkdir -p {params.out_dir}
        python workflow/scripts/aavf_filter.py \
            -i {input} \
            -o {output} \
            {params.flag} {params.value} \
            2> {log}
        """

rule aavf_filtered_vg:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq.aavf",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/lofreq_{{lofreq_value}}.aavf",
    params:
        out_dir=lambda wc: f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{wc.sample_id}/variants_filtered/{wc.lofreq_group}",
        flag=lambda wc: "--min-acf" if wc.lofreq_group == "AF" else "--min-acc",
        value=lambda wc: wc.lofreq_value
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/lofreq_{{lofreq_value}}.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants_filtered/{{lofreq_group}}/lofreq_{{lofreq_value}}.txt"
    shell:
        """
        mkdir -p {params.out_dir}
        python workflow/scripts/aavf_filter.py \
            -i {input} \
            -o {output} \
            {params.flag} {params.value} \
            2> {log}
        """

rule aavf_to_hivdb_linear:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/variants/alns_to_{{reference_name}}.aavf"
    output:
        f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/drug_resistance_report/alns_to_{{reference_name}}_results.txt"
    params:
        out_dir = f"{OUT_DIR}/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/drug_resistance_report"
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/drug_resistance_report/aavf_to_hivdb_alns_to_{{reference_name}}.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/{ALIGNER}/{{sample_id}}/drug_resistance_report/aavf_to_hivdb_alns_to_{{reference_name}}.txt"
    shell:
        """
        mkdir -p {params.out_dir} &&
        python workflow/scripts/aavf_to_hivdb.py --aavf {input} --out {output} 2> {log}
        """

rule aavf_to_hivdb_vg:
    input:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/variants/lofreq.aavf",
    output:
        f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/drug_resistance_report/results.txt",
    params:
        out_dir = f"{OUT_DIR}/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/drug_resistance_report",
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"logs/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/drug_resistance_report/aavf_to_hivdb_vg.log"
    benchmark:
        f"benchmarks/{BATCH_NAME}/vg/{ALIGNER}/{{sample_id}}/drug_resistance_report/aavf_to_hivdb_vg.txt"
    shell:
        """
        mkdir -p {params.out_dir} &&
        python workflow/scripts/aavf_to_hivdb.py --aavf {input} --out {output} 2> {log}
        """