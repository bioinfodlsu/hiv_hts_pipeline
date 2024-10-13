rule index_bowtie2:
    input:
       reference = config["reference"]
    output:
       touch("{0}".format(config["out_dir"])+"/bowtie2_index/index.done")
    params:
       index_basename="{0}".format(config["out_dir"])+"/bowtie2_index/index"
    conda:
       "../envs/bowtie2.yaml"
    shell:
       "bowtie2-build {input.reference} {params.index_basename}"

ruleorder: bowtie2_to_bam > align_bowtie2

rule align_bowtie2:
    input:
        reference_flag = "{0}/bowtie2_index/index.done".format(config["out_dir"]),
        sub1 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_1_subsampled.fq",
        sub2 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_2_subsampled.fq"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        index_basename = rules.index_bowtie2.params.index_basename,
        preset = lambda wildcards: config["aligner_params_dict"][wildcards.param_group]["preset"]
    conda:
       "../envs/bowtie2.yaml"
    shell:
        "bowtie2 --local --{params.preset} -x {params.index_basename}  -1 {input.sub1} -2 {input.sub2} | samtools view -bS - > {output}"

rule align_bowtie2_SE:
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/bowtie2_index/index.done",
        # query = config["reads"][list(config["reads"].keys())[0]][0]
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        index_basename = rules.index_bowtie2.params.index_basename,
        preset = lambda wildcards: config["aligner_params_dict"][wildcards.param_group]["preset"]
    conda:
       "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 --local --{params.preset} -x {params.index_basename}  -U {input.query} -S {output}
        """

rule bowtie2_to_bam:
    input:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.bam"
    threads:
        workflow.cores/len(config["reads"])
    params:
        index_basename = rules.index_bowtie2.params.index_basename,
        preset = lambda wildcards: config["aligner_params_dict"][wildcards.param_group]["preset"]
    conda:
       "../envs/bowtie2.yaml"
    shell:
        """
        picard SamFormatConverter I={input} O={output}
        """