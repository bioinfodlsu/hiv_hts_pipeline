#from snakemake.utils import min_version
#min_version("6.0")
#
# module last_workflow:
#     snakefile:"rules/last.smk"
#     config:config
# use rule * from last_workflow
#
# module bowtie_workflow:
#     snakefile:"rules/bowtie.smk"
#     config:config
# use rule* from bowtie_workflow
from snakemake.utils import Paramspace
import pandas as pd

last_params_pd = pd.read_csv(config["last_param_space_file"],sep="\t")
last_paramspace = Paramspace(last_params_pd)
bowtie2_params_pd = pd.read_csv(config["bowtie2_param_space_file"],sep="\t")
bowtie2_paramspace = Paramspace(bowtie2_params_pd)

include:"rules/last.smk"
include:"rules/bowtie.smk"

rule all: #Initial Rule in snakemake
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns",
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns",
        "{0}".format(config['out_dir'])+"/last_trained_alignments/last_counts_alns"

        #expand("{out_dir}/bowtie2_alignments/{sample_id}.{params}/alns.sam",
        #         out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=bowtie2_paramspace.instance_patterns)
        # expand("{out_dir}/last_alignments/{sample_id}/{params}/alns.tab",
        #          out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=last_paramspace.instance_patterns)
        #expand("{out_dir}/last_alignments/{sample_id}.tab",
        #         out_dir = config["out_dir"],sample_id=config["reads"].keys())
rule count_last_trained_alns:
    input:
        aln_files = expand("{out_dir}/last_trained_alignments/{sample_id}.tab",
                           out_dir = config["out_dir"],sample_id=config["reads"].keys())
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/last_counts_alns"
    run:
        out = open(output[0],"w")
        out.write("samples\t"+"\taln_counts\n")
        for aln_filename in input.aln_files:
            with open(aln_filename) as f:
                aln_count = 0
                for line in f:
                    if not line.startswith("#"):
                        aln_count +=1

            #separate aln_filename into the wildcards and params
            filename_split = aln_filename.split("/")
            index = filename_split.index("last_trained_alignments")
            sample_id = filename_split[index+1]
            towrite = sample_id+"\t"+f"{aln_count}"+"\n"
            out.write(towrite)
        out.close()

rule count_last_alns:
    input:
        aln_files=expand("{out_dir}/last_alignments/{sample_id}/{params}/alns.tab",
                         out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=last_paramspace.instance_patterns)
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns"
    run:
        out = open(output[0],"w")
        param_names = list(last_params_pd)
        out.write("samples\t"+".".join(param_names)+"\taln_counts\n")
        for aln_filename in input.aln_files:
            with open(aln_filename) as f:
                aln_count = 0
                for line in f:
                    if not line.startswith("#"):
                        aln_count +=1

            #separate aln_filename into the wildcards and params
            filename_split = aln_filename.split("/")
            index = filename_split.index("last_alignments")
            sample_id = filename_split[index+1]
            params_value=filename_split[index+2:-1]
            values = ".".join([ pair.split("~")[-1] for pair in params_value])
            towrite = sample_id+"\t"+values+"\t"+f"{aln_count}"+"\n"
            out.write(towrite)
        out.close()



rule count_bowtie2_alns:
    input:
        aln_files= expand("{out_dir}/bowtie2_alignments/{sample_id}/{params}/alns.sam",
                           out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=bowtie2_paramspace.instance_patterns)
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns"
    run:
        out = open(output[0],"w")
        param_names = list(bowtie2_params_pd)
        out.write("samples\t"+".".join(param_names)+"\taln_counts\n")
        for aln_filename in input.aln_files:
            with open(aln_filename) as f:
                aln_count = 0
                for line in f:
                    if not line.startswith("#"):
                        aln_count +=1

            #separate aln_filename into the wildcards and params
            filename_split = aln_filename.split("/")
            index = filename_split.index("bowtie2_alignments")
            sample_id = filename_split[index+1]
            params_value=filename_split[index+2:-1]
            values = ".".join([ pair.split("~")[-1] for pair in params_value])
            towrite = sample_id+"\t"+values+"\t"+f"{aln_count}"+"\n"
            out.write(towrite)
        out.close()
