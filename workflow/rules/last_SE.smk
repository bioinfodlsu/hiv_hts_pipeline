from os import path
import itertools

rule last_frag_stat_est_single: #Rule for estimating length distribution of paired-end reads
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index_single.done",
        # query = config["reads"][list(config["reads"].keys())[0]][0]
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq"
    params:
        last_index_basename= "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index",
        last_al = "{0}".format(config["out_dir"])+"/last_frag_stat".format(config["out_dir"])
    output:
        frag_len_est="{0}".format(config["out_dir"])+"/last_frag_stat/{sample_id}_single.frag_len_est"
    conda:
       "../envs/aavf.yaml"
    shell:
        #note that the lastal for last-pair-probs uses LAST's default scoring scheme instead of learned scores
        """
        set +o pipefail;
        mkfifo {params.last_al}/{wildcards.sample_id} |
        gzip -cdf {input.query} | head -n 400000 > {params.last_al}/{wildcards.sample_id} &
        fastq-interleave {params.last_al}/{wildcards.sample_id} | lastal -Q1 -i1 {params.last_index_basename} > {output.frag_len_est}
        rm {params.last_al}/{wildcards.sample_id}
        """

rule align_last_trained_SE:
# Rule for aligning paired-end reads to a reference genome, with score parameters provided by user.
# It would have been nicer to incorporate this in the rule align_last, but it is too messy and hurts readability and debuggability
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        scoring = "{0}".format(config["out_dir"])+"/last_score_sample/{sample_id}_scoring_scheme",
        # probably using the wrong query sequence here
        # TODO: change query to the correct query sequence
        # query = config["reads"][list(config["reads"].keys())[0]][0],
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq",
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        # mafsam = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.sam"
    params:
        last_index_basename="{0}".format(config["out_dir"])+"/last_index/{sample_id}_index",
        mean = getMean,
        std = getSTD
    threads:
        workflow.cores/len(config["reads"])
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.maf"
    conda:
       "../envs/last.yaml"
    shell:
        """
        lastal -Q1 -i1 -p {input.scoring} {params.last_index_basename} {input.query} |
        last-map-probs -m 0.95 > {output}
        """

rule align_last_trained_SE_toSAM:
# Rule for aligning paired-end reads to a reference genome, with score parameters provided by user.
# It would have been nicer to incorporate this in the rule align_last, but it is too messy and hurts readability and debuggability
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        scoring = "{0}".format(config["out_dir"])+"/last_score_sample/{sample_id}_scoring_scheme",
        # probably using the wrong query sequence here
        # TODO: change query to the correct query sequence
        # query = config["reads"][list(config["reads"].keys())[0]][0],
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq",
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        maf = "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.maf"
    params:
        last_index_basename="{0}".format(config["out_dir"])+"/last_index/{sample_id}_index",
        mean = getMean,
        std = getSTD
    threads:
        workflow.cores/len(config["reads"])
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    conda:
       "../envs/last.yaml"
    shell:
        """
        maf-convert sam {input.maf} -f {input.dict} > {output}
        """
    
rule align_last_SE: #Rule for aligning paired-end reads to a reference genome, with score training
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        # query = config["reads"][list(config["reads"].keys())[0]][0],
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq",
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        # mafsam = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.sam"
    params:
        last_index_basename="{0}".format(config["out_dir"])+"/last_index/{sample_id}_index",
        mean = getMean,
        std = getSTD,
        fromParamsFile = lambda wildcards: config["last_params_dict"][wildcards.param_group]
        # fromParamsFile=last_paramspace.instance
    threads:
        workflow.cores/len(config["reads"])
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.maf"
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        lastal -Q1 -i1 -r {params.fromParamsFile[match]} -q {params.fromParamsFile[mismatch]} -a {params.fromParamsFile[gapOpen]} -b {params.fromParamsFile[gapExtend]} {params.last_index_basename} {input.query} |
        last-map-probs -m 0.95 > {output}
        """

rule align_last_SE_toSAM: #Rule for aligning paired-end reads to a reference genome, with score training
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        # query = config["reads"][list(config["reads"].keys())[0]][0],
        query = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_subsampled.fq",
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        maf = "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.maf"
    params:
        last_index_basename="{0}".format(config["out_dir"])+"/last_index/{sample_id}_index",
        mean = getMean,
        std = getSTD,
        fromParamsFile = lambda wildcards: config["last_params_dict"][wildcards.param_group]
        # fromParamsFile=last_paramspace.instance
    threads:
        workflow.cores/len(config["reads"])
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        maf-convert sam {input.maf} -f {input.dict} > {output}
        """
