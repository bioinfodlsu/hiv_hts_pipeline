rule convertSAM:
    input:
        reference = config["reference"],
        sam = "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        last_index_basename = "{0}/last_index/index".format(config["out_dir"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        samtools faidx {input.reference}
        samtools view -bt {params.last_index_basename} {input.sam} > {output}
        """

rule convertSAM_trained:
    input:
        reference = config["reference"],
        sam = "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        last_index_basename = "{0}/last_index/index".format(config["out_dir"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        samtools faidx {input.reference}
        samtools view -bt {params.last_index_basename} {input.sam} > {output}
        """

rule convertBAM:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.bam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/blast.yaml"
    shell:
        """
        samtools sort -n {input.bam} |
        samtools view -h - > {output}
        """

rule sortBAM:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.bam"
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    conda:
        "../envs/last.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule sortSAM:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sorted.sam"
    conda:
        "../envs/last.yaml"
    shell:
        "samtools sort -n {input} > {output}"

rule sortSAM_trained:
    input:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sorted.sam"
    conda:
        "../envs/last.yaml"
    shell:
        "samtools sort -n {input} > {output}"

rule indexBAM:
    input:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sorted.bam"
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sorted.bam.bai"
    conda:
        "../envs/{aligner}.yaml"
    shell:
        "samtools index {input}"
