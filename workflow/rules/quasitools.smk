rule create_bed:
    input:
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        reference = config["reference"],
        annotation = "test_data/hxb2_gene_annotation.bed"
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}.bed"
    params:
        original_bed = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}_og.bed",
        annotated = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}_annotated.bed"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        bedtools bamtobed -i {input.bam} > {params.original_bed} &&
        bedtools intersect -a {params.original_bed} -b {input.annotation} -wa -wb > {params.annotated} &&
        awk '{{print $1 "\t" $2 "\t" $3 "\t" $10}}' {params.annotated} > {output}
        """

rule aavf_quasitools:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bai",
        # bed = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}.bed"
        bed = "test_data/hxb2_gene_annotation.bed",
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_filled.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_quasitools.aavf"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        quasitools call aavar {input.bam} {input.reference} {input.bed} -o {output}
        """
