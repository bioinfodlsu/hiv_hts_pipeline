rule variants_lofreq:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bai"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        lofreq faidx {input.reference}
        lofreq call-parallel --pp-threads {threads} -f {input.reference} -o {output} {input.bam}
        """

rule bcftools_index:
    input:
        reference = config["reference"],
        dictionary = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_partial.vcf",
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        picard UpdateVcfSequenceDictionary -I {input.vcf} -O {output} -SD {input.dictionary}
        """

rule bcftools_fillAC:
    input:
        reference = config["reference"],
        dictionary = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_partial.vcf",
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_filled.vcf"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ob -o {output} -- -t AC
        """

# rule variants_lofreq_filtered:
#     input:
#         reference = config["reference"],
#         bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
#         bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bai"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.vcf"
#     params:
#         e = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         """
#         lofreq faidx {input.reference}
#         lofreq call-parallel -q {params.e[QUAL]} --pp-threads {threads} -f {input.reference} -o {output} {input.bam}
#         """

