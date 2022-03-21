rule convertSAM:
    input:
        reference = config["reference"],
        sam = "{0}".format(config['out_dir'])+"/last_alignments/{sample_id}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{sample_id}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        last_index_basename = "{0}/last_index/index".format(config["out_dir"])
    conda:
        "../envs/{aligner}.yaml"
    shell:
        """
        samtools faidx {input.reference}
        samtools view -bt {params.last_index_basename} {input.sam} > {output}
        """

rule sortBAM:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{sample_id}/paramgroup_{param_group}/alns.bam"
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    conda:
        "../envs/{aligner}.yaml"
    shell:
        "samtools sort -o {output} {input}"

rule indexBAM:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam.bai"
    conda:
        "../envs/{aligner}.yaml"
    shell:
        "samtools index {input}"
