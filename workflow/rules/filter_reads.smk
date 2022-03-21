rule trim_setqk:
    input:
        reads1 = lambda wildcards: config["reads"][wildcards.sample_id][0],
        reads2 = lambda wildcards: config["reads"][wildcards.sample_id][1]
    output:
        trim1 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_1_trimmed.fq",
        trim2 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_2_trimmed.fq"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/{}.yaml".format(config['aligner'])
    shell:
        """
        seqtk trimfq {input.reads1} > {output.trim1}
        seqtk trimfq {input.reads2} > {output.trim2}
        """

rule subsample_seqtk:
    input:
        trim1 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_1_trimmed.fq",
        trim2 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_2_trimmed.fq"
    output:
        sub1 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_1_subsampled.fq",
        sub2 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}_2_subsampled.fq"
    params:
        n = 1000000,      # number of reads after subsampling
        seed = 100      # seed to initialize a pseudorandom number generator
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/{}.yaml".format(config['aligner'])
    shell:
        """
        seqtk sample -s {params.seed} {input.trim1} {params.n} > {output.sub1}
        seqtk sample -s {params.seed} {input.trim2} {params.n} > {output.sub2}
        """
