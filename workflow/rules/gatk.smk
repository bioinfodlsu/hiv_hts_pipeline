# rule create_dict_gatk:
#     input:
#         reference = config["reference"]
#     output:
#        "test_data/hxb2.dict"
#     conda:
#        "../envs/gatk.yaml"
#     shell:
#        "gatk CreateSequenceDictionary -R {input.reference} -O {output}"

rule add_replace_read_groups:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.merged.bam",
    params:
        rgid_id = "1",
        rglb_string = "lib1",
        rgpl_string = "illumina",
        rgsm_string = "sample1",
        rgpu_string = "illumina"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        """
            gatk AddOrReplaceReadGroups -I {input} \
            -O {output} \
            -RGID {params.rgid_id} \
            -RGLB {params.rglb_string} \
            -RGPL {params.rgpl_string} \
            -RGSM {params.rgsm_string} \
            -RGPU {params.rgpu_string}
        """

rule sort_gatk:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.merged.bam",
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.merged.sorted.bam",
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SortSam -I {input} -O {output} -SO coordinate"


rule mark_duplicates:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.merged.sorted.bam",
    output:
        md = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.marked_duplicates.bam",
        mdm = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.metrics"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk MarkDuplicates I={input.bam} O={output.md} M={output.mdm}"

rule split_n_trim:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.marked_duplicates.bam",
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.split.bam",
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk SplitNCigarReads -R {input.reference} -I {input.bam} -O {output}"

rule bam_index_gatk:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.split.bam",
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.split.bam.bai",
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        "samtools index {input} > {output}"

rule variants_gatk:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.split.bam",
        bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.split.bam.bai",
    output:
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk HaplotypeCaller -R {input.reference} -I {input.bam} -O {output}"

rule variants_gatk_filtered:
    input:
        reference = config["reference"],
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/filtered/variants_{aligner}_{lofreq_param_group}.vcf"
    params:
        e = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk VariantFiltration  -V {input.vcf} \
        -filter "QUAL < {params.e[QUAL]}" \
        --filter-name "QUAL{params.e[QUAL]}" \
        -O {output}
        """