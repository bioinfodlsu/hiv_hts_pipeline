
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
    threads:
        workflow.cores/len(config["reads"])
    params:
        index_basename = rules.index_bowtie2.params.index_basename
    output:
        "{0}/bowtie2_alignments/{{sample_id}}.sam".format(config["out_dir"])
    shell:
        "bowtie2 --local -x {params.index_basename}  -1 {input.reads1} -2 {input.reads2} -S {output}"
