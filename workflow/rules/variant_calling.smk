rule variants_lofreq:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam.bai"
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

# add rule to convery VCF to AAVF 