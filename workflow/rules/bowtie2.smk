
rule index_bowtie2:
    input:
       reference = config["reference"]
    output:
       touch("{0}/bowtie2_index/index.done".format(config["out_dir"]))
    params:
       index_basename="{0}/bowtie2_index/index".format(config["out_dir"])
    shell:
       "bowtie2-build {input.reference} {params.index_basename}"

rule align_bowtie2:
    input:
        reference_flag = "{0}/bowtie2_index/index.done".format(config["out_dir"]),
        reads1 = lambda wildcards: config["reads"][wildcards.sample_id][0],
        reads2 = lambda wildcards: config["reads"][wildcards.sample_id][1]
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        index_basename = rules.index_bowtie2.params.index_basename,
        preset = lambda wildcards: config["bowtie2_params_dict"][wildcards.param_group]["preset"],
        X = lambda wildcards: config["bowtie2_params_dict"][wildcards.param_group]["X"]

    shell:
        #"{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.sam"
        #"{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.sam | samtools view -bS - > {output}"
        "bowtie2 --local --{params.preset} -x {params.index_basename}  -1 {input.reads1} -2 {input.reads2} | samtools view -bS - > {output}"
