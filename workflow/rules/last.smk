from os import path
import itertools

if "filetype" in config:
    if config["filetype"] == "fasta":
        interleave = srcdir("scripts/fasta-interleave.sh")
        file_type = "fasta"
    elif config["filetype"] == "fastq":
        interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
        file_type = "fastq"
else:
    interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
    file_type = "fastq"

#Input functions
def get_read1(wildcards):
    return "{}/filtered_reads/{}/{}_1_subsampled.fq".format(config['out_dir'],wildcards.sample_id,wildcards.sample_id);

def get_read2(wildcards):
    return "{}/filtered_reads/{}/{}_2_subsampled.fq".format(config['out_dir'],wildcards.sample_id,wildcards.sample_id);

def getMean(wildcards):
    est_path = "{}/last_frag_stat/{}.frag_len_est".format(config["out_dir"],wildcards.sample_id)
    if(path.exists(est_path)):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated mean distance"):
                    return(line.rstrip().split()[-1])
    else:
        return(-1)

def getSTD(wildcards):
    est_path = "{}/last_frag_stat/{}.frag_len_est".format(config["out_dir"],wildcards.sample_id)
    if(path.exists(est_path)):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated standard deviation of distance"):
                    return(line.rstrip().split()[-1])
    else:
        return(-1)

def get_input_align_last(wildcards):
    input_dict = {}
    input_dict["reference_flag"] = "{0}/last_index/index.done".format(config["out_dir"])
    input_dict["reads1"] = get_read1(wildcards)
    input_dict["reads2"] = get_read2(wildcards)
    return input_dict


#Rule declarations begin here:
rule last_all: #Initial Rule in snakemake
    input:
        "{0}/last_score_sample/scoring_scheme".format(config["out_dir"])
        #expand("{out_dir}/last_alignments/{sample_id}.tab", out_dir = config["out_dir"],sample_id=config["reads"].keys())

rule last_db: #Rule for constructing LAST index
    input:
        reference = config["reference"]
    output:
       touch("{0}/last_index/index.done".format(config["out_dir"]))
    params:
       index_basename="{0}/last_index/index".format(config["out_dir"])
    conda:
       "../envs/last.yaml"
    shell:
       "lastdb {params.index_basename} {input.reference}"

rule last_score_training: #Rule for training score parameters
    input:
        reference_flag = "{0}/last_index/index.done".format(config["out_dir"]),
        query = config["reads"][list(config["reads"].keys())[0]][0]
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
        last_sample1="{0}/last_score_sample/sample_1.sample.fa".format(config["out_dir"]),
        file_type=file_type
    output:
        "{0}/last_score_sample/scoring_scheme".format(config["out_dir"])
    conda:
       "../envs/last.yaml"
    shell:
        """
        set +o pipefail;
        gzip -cdf {input.query} | head -n 400000 > {params.last_sample1}
        last-train --verbose --revsym --matsym --gapsym --sample-number=400000 -Q1  {params.last_index_basename} {params.last_sample1} > {output}
        """

rule last_frag_stat_est: #Rule for estimating length distribution of paired-end reads
    input:
        unpack(get_input_align_last)
    params:
        last_index_basename= "{0}/last_index/index".format(config["out_dir"]),
        last_al = "{0}/last_frag_stat".format(config["out_dir"])
    output:
        frag_len_est="{0}".format(config["out_dir"])+"/last_frag_stat/{sample_id}.frag_len_est"
    conda:
       "../envs/last.yaml"
    shell:
        #note that the lastal for last-pair-probs uses LAST's default scoring scheme instead of learned scores
        """
        set +o pipefail;
        mkfifo {params.last_al}/{wildcards.sample_id}_1 |
        mkfifo {params.last_al}/{wildcards.sample_id}_2 |
        gzip -cdf {input.reads1} | head -n 400000 > {params.last_al}/{wildcards.sample_id}_1 &
        gzip -cdf {input.reads2} | head -n 400000 > {params.last_al}/{wildcards.sample_id}_2 &
        fastq-interleave {params.last_al}/{wildcards.sample_id}_1 {params.last_al}/{wildcards.sample_id}_2 | lastal -Q1 -i1 {params.last_index_basename} | last-pair-probs -e > {output.frag_len_est}
        rm {params.last_al}/{wildcards.sample_id}_1 {params.last_al}/{wildcards.sample_id}_2
        """

rule create_dict:
    input:
        reference = config["reference"]
    output:
       "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict"
    conda:
       "../envs/last.yaml"
    shell:
       "picard CreateSequenceDictionary R={input.reference} O={output}"

rule align_last_trained:
# Rule for aligning paired-end reads to a reference genome, with score parameters provided by user.
# It would have been nicer to incorporate this in the rule align_last, but it is too messy and hurts readability and debuggability
    input:
        unpack(get_input_align_last),
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        scoring = "{0}/last_score_sample/scoring_scheme".format(config["out_dir"]),
        frag_len_est = "{0}".format(config["out_dir"])+"/last_frag_stat/{sample_id}.frag_len_est",
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
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
        fastq-interleave {input.reads1} {input.reads2} |
        lastal -Q1 -i1 -p {input.scoring} {params.last_index_basename} |
        last-pair-probs -f {params.mean} -s {params.std} -m 0.01 -d 0.1 |
        maf-convert sam -f {input.dict} > {output}
        """
    
rule align_last: #Rule for aligning paired-end reads to a reference genome, with score training
    input:
        unpack(get_input_align_last),
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict",
        frag_len_est = "{0}".format(config["out_dir"])+"/last_frag_stat/{sample_id}.frag_len_est"
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
        mean = getMean,
        std = getSTD,
        fromParamsFile = lambda wildcards: config["last_params_dict"][wildcards.param_group]
        # fromParamsFile=last_paramspace.instance
    threads:
        workflow.cores/len(config["reads"])
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam"
    conda:
        "../envs/last.yaml"
    shell:
        """
        fastq-interleave {input.reads1} {input.reads2} |
        lastal -Q1 -i1 -r {params.fromParamsFile[match]} -q {params.fromParamsFile[mismatch]} -a {params.fromParamsFile[gapOpen]} -b {params.fromParamsFile[gapExtend]} {params.last_index_basename} |
        last-pair-probs -f {params.mean} -s {params.std} -m 0.01 -d 0.1 |
        maf-convert sam -f {input.dict} > {output}
        """