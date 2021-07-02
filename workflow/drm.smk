#Pipeline for DRM calling
#Here we will read config file which contains information about which tools to use
#for , and write rules to decide will rules in the rules folder to call.
from snakemake.utils import min_version
min_version("6.0")

from snakemake.utils import Paramspace
import pandas as pd
bowtie2_params_pd = pd.read_csv(config["bowtie2_param_space_file"],sep="\t")
bowtie2_paramspace = Paramspace(bowtie2_params_pd)

module module_bowtie2:
    snakefile:"rules/bowtie.smk"
    config:config

use rule * from module_bowtie2


rule all:
    input:
        expand("{out_dir}/bowtie2_alignments/{sample_id}/{params}/alns.sam",
                           out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=bowtie2_paramspace.instance_patterns)
