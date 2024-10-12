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
        lofreq call-parallel --no-default-filter -N --pp-threads {threads} -f {input.reference} -o {output} {input.bam}
        """
        # lofreq call-parallel --pp-threads {threads} -f {input.reference} -o {output} {input.bam}

rule variants_lofreq_filtered:
    input:
        reference = config["reference"],
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants_filtered/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.vcf"
    params:
        e = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/lofreq.yaml"
    shell:
        # """
        # lofreq filter -i {input.vcf} --verbose -Q {params.e[QUAL]} -o {output}
        # """
        """
        lofreq filter -i {input.vcf} --verbose -a {params.e[AF]} -o {output}
        """
        # """
        # lofreq filter -i {input.vcf} --verbose -V {params.e[DP]} -o {output}
        # """