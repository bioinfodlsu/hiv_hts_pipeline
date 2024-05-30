from os import path
import itertools

# if "filetype" in config:
#     if config["filetype"] == "fasta":
#         interleave = srcdir("scripts/fasta-interleave.sh")
#         file_type = "fasta"
#     elif config["filetype"] == "fastq":
#         interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
#         file_type = "fastq"
# else:
#     interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
#     file_type = "fastq"

#Input functions
# def get_read(wildcards):
#     return "{}/filtered_reads/{}/{}_1_subsampled.fq".format(config['out_dir'],wildcards.sample_id,wildcards.sample_id);

# def getMean(wildcards):
#     est_path = "{}/last_frag_stat/{}.frag_len_est".format(config["out_dir"],wildcards.sample_id)
#     if(path.exists(est_path)):
#         with open(est_path) as f:
#             for line in f:
#                 if line.startswith("# estimated mean distance"):
#                     return(line.rstrip().split()[-1])
#     else:
#         return(-1)

# def getSTD(wildcards):
#     est_path = "{}/last_frag_stat/{}.frag_len_est".format(config["out_dir"],wildcards.sample_id)
#     if(path.exists(est_path)):
#         with open(est_path) as f:
#             for line in f:
#                 if line.startswith("# estimated standard deviation of distance"):
#                     return(line.rstrip().split()[-1])
#     else:
#         return(-1)

# def get_input_align_last(wildcards):
#     input_dict = {}
#     input_dict["reference_flag"] = "{0}/last_index/index.done".format(config["out_dir"])
#     input_dict["reads"] = get_read(wildcards)
#     return input_dict


rule last_frag_stat_est_single: #Rule for estimating length distribution of paired-end reads
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index_single.done",
        query = config["reads"][list(config["reads"].keys())[0]][0]
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

# rule create_dict:
#     input:
#         reference = config["reference"]
#     output:
#        "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict"
#     conda:
#        "../envs/last.yaml"
#     shell:
#        "picard CreateSequenceDictionary R={input.reference} O={output}"
# ruleorder: align_last_trained > align_last_trained_single

rule align_last_trained_single:
# Rule for aligning paired-end reads to a reference genome, with score parameters provided by user.
# It would have been nicer to incorporate this in the rule align_last, but it is too messy and hurts readability and debuggability
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        scoring = "{0}".format(config["out_dir"])+"/last_score_sample/{sample_id}_scoring_scheme",
        query = config["reads"][list(config["reads"].keys())[0]][0],
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        # mafsam = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.sam"
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
        lastal -Q1 -i1 -p {input.scoring} {params.last_index_basename} {input.query} |
        last-map-probs -m 0.95 |
        maf-convert sam -f {input.dict} > {output}
        """
    
rule align_last_single: #Rule for aligning paired-end reads to a reference genome, with score training
    input:
        reference_flag = "{0}".format(config["out_dir"])+"/last_index/{sample_id}_index.done",
        query = config["reads"][list(config["reads"].keys())[0]][0],
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
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        lastal -Q1 -i1 -r {params.fromParamsFile[match]} -q {params.fromParamsFile[mismatch]} -a {params.fromParamsFile[gapOpen]} -b {params.fromParamsFile[gapExtend]} {params.last_index_basename} {input.query} |
        last-map-probs -m 0.95 |
        maf-convert sam -f {input.dict} > {output}
        """

rule create_bed:
    input:
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        reference = config["reference"],
        annotation = "test_data/hxb2_gene_annotation.bed"
    output:
        "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}.bed"
    params:
        original_bed = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}_og.bed",
        annotated = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}_annotated.bed"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        bedtools bamtobed -i {input.bam} > {params.original_bed} &&
        bedtools intersect -a {params.original_bed} -b {input.annotation} -wa -wb > {params.annotated} &&
        awk '{{print $1 "\t" $2 "\t" $3 "\t" $10}}' {params.annotated} > {output}
        """

rule aavf_quasitools:
    input:
        reference = config["reference"],
        bam = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bam",
        bai = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sorted.bai",
        # bed = "{0}".format(config['out_dir'])+"/{aligner}_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/{sample_id}.bed"
        bed = "test_data/hxb2_gene_annotation.bed",
        vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_filled.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_quasitools.aavf"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/aavf.yaml"
    shell:
        """
        quasitools call aavar {input.bam} {input.reference} {input.bed} -o {output}
        """
