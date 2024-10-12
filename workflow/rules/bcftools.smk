# rule bcftools_index:
#     input:
#         reference = config["reference"],
#         dictionary = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
#         bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_partial.vcf",
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         """
#         picard UpdateVcfSequenceDictionary -I {input.vcf} -O {output} -SD {input.dictionary}
#         """

# rule bcftools_fillAC:
#     input:
#         reference = config["reference"],
#         dictionary = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
#         bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_partial.vcf",
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_filled.vcf"
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         """
#         bcftools +fill-tags {input.vcf} -Ob -o {output} -- -t AC
#         """

# rule variants_bcftools:
#     input:
#         reference = config["reference"],
#         bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
#         bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bai"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         """
#         bcftools mpileup -Ou -f {input.reference} {input.bam} | bcftools call -mv -Ob -o {output}
#         """

# rule variants_bcftools_filtered:
#     input:
#         reference = config["reference"],
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
#     output:
#         "{0}".format(config['out_dir'])+"/variants_filtered/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.vcf"
#     params:
#         e = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         # """
#         # lofreq filter -i {input.vcf} --verbose -Q {params.e[QUAL]} -o {output}
#         # """
#         # """
#         # lofreq filter -i {input.vcf} --verbose -A {params.e[AF]} -o {output}
#         # """
#         """
#         lofreq filter -i {input.vcf} --verbose -V {params.e[DP]} -o {output}
#         """

# rule bcftools_fillAF:
#     input:
#         reference = config["reference"],
#         dictionary = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
#         bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_partial.vcf",
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_filled.vcf"
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         """
#         bcftools +fill-tags {input.vcf} -Ob -o {output} -- -t AF
#         """