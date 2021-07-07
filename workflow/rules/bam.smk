rule sortBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"


rule indexBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam.bai"
    shell:
        "samtools index {input}"
