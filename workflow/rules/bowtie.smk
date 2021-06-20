
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
        fromFile=bowtie2_paramspace.instance,
        index_basename = rules.index_bowtie2.params.index_basename,

    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}"+f".{bowtie2_paramspace.wildcard_pattern}.sam"
        #Snakemake is confusing. note that in each of the components above, the curly braces serve a different purpose.
    shell:
        "bowtie2 --local --{params.fromFile[preset]} -x {params.index_basename}  -1 {input.reads1} -2 {input.reads2} -S {output}"
        #Snakemake is confusing. Why is preset allowed but not 'preset' ??
