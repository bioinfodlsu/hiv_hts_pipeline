# Do we need this? Maybe directly sort?

rule sortBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.sorted.bam"
    shell:
        samtools sort -o {output} {input}


rule indexBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.sorted.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.sorted.bam.bai"
    shell:
        samtools index {input}
