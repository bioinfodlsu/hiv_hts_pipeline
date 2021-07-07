rule lofreq:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam",
        bai = "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam.bai"

    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/variants.vcf"

    threads: 10

    shell:
        """
        lofreq faidx {input.reference}
        lofreq call-parallel --pp-threads {threads} -f {input.reference} -o {output} {input.bam}
        """
