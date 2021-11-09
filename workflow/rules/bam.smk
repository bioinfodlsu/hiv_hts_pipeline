rule sortBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "samtools sort -o {output} {input}"


rule indexBAM:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam.bai"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "samtools index {input}"
