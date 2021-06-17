from os import path
import itertools
#configfile: "testdata/config_test.yaml"

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
    return config["reads"][wildcards.sample_id][0]

def get_read2(wildcards):
    return config["reads"][wildcards.sample_id][1]


def getMean(wildcard):
    est_path = "{}/last_alignments/{}.frag_len_est".format(config["out_dir"],wildcard)
    if(path.exists(est_path)):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated mean distance"):
                    return(line.rstrip().split()[-1])
    else:
        return(-1)

def getSTD(wildcard):
    est_path = "{}/last_alignments/{}.frag_len_est".format(config["out_dir"],wildcard)
    if(path.exists(est_path)):
        with open(est_path) as f:
            for line in f:
                if line.startswith("# estimated standard deviation of distance"):
                    return(line.rstrip().split()[-1])
    else:
        return(-1)

def getScoring(wildcards):
    if "training" in config:
        if config["training"].lower() != "no":
            return "{}/score_sample/scoring_scheme".format(config["out_dir"])
        else: return "HOXD70"
    else:
        return "HOXD70"

def get_input_align_last(wildcards):
    input_dict = {}
    input_dict["reference_flag"] = "{}/last_index/index.done".format(config["out_dir"])
    input_dict["reads1"] = get_read1(wildcards)
    input_dict["reads2"] = get_read2(wildcards)
    #input_dict["frag_len_est"] = config["out_dir"]+"last_alignments/{sample_id}.frag_len_est"
    if "training" in config:
        if config["training"].lower() != "no":
            input_dict["scoring"] = "{0}/score_sample/scoring_scheme".format(config["out_dir"])
    else:
        input_dict["scoring"] = "{0}/score_sample/scoring_scheme".format(config["out_dir"])
    return input_dict


# Rule declarations begin here:
rule last_all: #Initial Rule in snakemake
    input:
        "{0}/score_sample/scoring_scheme".format(config["out_dir"])
        #expand("{out_dir}/last_alignments/{sample_id}.tab", out_dir = config["out_dir"],sample_id=config["reads"].keys())

rule last_db: #Rule for constructing LAST index
    input:
        reference = config["reference"]
    output:
       touch("{0}/last_index/index.done".format(config["out_dir"]))
    params:
       index_basename="{0}/last_index/index".format(config["out_dir"])
    #conda:
    #    "env/last.yaml"
    shell:
       "lastdb {params.index_basename} {input.reference}"

rule last_score_training: #Rule for training score parameters
    input:
        reference_flag = "{0}/last_index/index.done".format(config["out_dir"]),
        query = config["reads"][list(config["reads"].keys())[0]][0]
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
        last_sample1="{0}/score_sample/sample_1.sample.fa".format(config["out_dir"]),
        file_type=file_type
    output:
        "{0}/score_sample/scoring_scheme".format(config["out_dir"])
    #conda:
    #    "env/lastxpy.yaml"
    shell:
        """
        set +o pipefail;
        gzip -cdf {input.query} | head -n 400000 > {params.last_sample1}
        last-train --verbose --revsym --matsym --gapsym --sample-number=400000 -Q1  {params.last_index_basename} {params.last_sample1} > {output}
        """
rule last_frag_stat_est: #Rule for estimating length distribution of paired-end reads
    input:
        unpack(get_input_align_last)
        #reference_flag = config["out_dir"]+"last_index/index.done",
        #reads1 = get_read1,
        #reads2 = get_read2,
        #scoring = config["out_dir"]+"last_sample/scoring_scheme",
    params:
        last_index_basename= "{0}/last_index/index".format(config["out_dir"]),
        last_al = "{0}/last_alignments".format(config["out_dir"]),
        scoring = getScoring
    output :
        frag_len_est="{0}/last_alignments/{{sample_id}}.frag_len_est".format(config["out_dir"])
    #conda:
    #    "env/last.yaml"
    shell:
        """
        set +o pipefail;
        mkfifo {params.last_al}/{wildcards.sample_id}_1
        mkfifo {params.last_al}/{wildcards.sample_id}_2
        gzip -cdf {input.reads1} | head -n 400000 > {params.last_al}/{wildcards.sample_id}_1 &
        gzip -cdf {input.reads2} | head -n 400000 > {params.last_al}/{wildcards.sample_id}_2 &
        fastq-interleave {params.last_al}/{wildcards.sample_id}_1 {params.last_al}/{wildcards.sample_id}_2 | lastal -Q1 -i1 -p {params.scoring}  {params.last_index_basename} | last-pair-probs -e > {output.frag_len_est}
        rm {params.last_al}/{wildcards.sample_id}_1 {params.last_al}/{wildcards.sample_id}_2
        """

rule align_last: #Rule for aligning RNA-seq reads to protein database
    input:
        unpack(get_input_align_last),
        frag_len_est = "{0}/last_alignments/{{sample_id}}.frag_len_est".format(config["out_dir"])
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
        mean = getMean,
        std = getSTD,
        #scoring = config["out_dir"]+"last_sample/scoring_scheme"
        scoring = getScoring
    threads: workflow.cores/len(config["reads"])
    output:
        "{0}/last_alignments/{{sample_id}}.tab".format(config["out_dir"])
    #conda:
    #    "env/lastal.yaml"
    shell:
        "fastq-interleave {input.reads1} {input.reads2} | "
        "parallel --gnu --pipe -L4 -j {threads} 'lastal -Q1 -i1 -p {params.scoring} {params.last_index_basename}' |"
        "last-pair-probs -f {params.mean} -s {params.std} -m 0.01 -d 0.1 |"
        "maf-convert tab > {output}"
