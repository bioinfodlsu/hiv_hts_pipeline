
# Parse Bowtie2 parameters.
#TODO Error checks
bowtie2_params = config["bowtie2_params"]
bowtie2_params_group_names = bowtie2_params[1:]
bowtie2_params_header = bowtie2_params[0]
bowtie2_params_dict = {}
for param_list in bowtie2_params[1:]:
    group_name = param_list[0]
    bowtie2_params_dict[group_name] = dict(zip(bowtie2_params_header,param_list))

config['bowtie2_params_dict'] = bowtie2_params_dict
#print(config['bowtie2_params_dict'] )


include: "rules/bowtie.params.smk"

rule all:
    input:
        expand("{out_dir}/bowtie2_alignments/{sample_id}/paramgroup_{param_group}.bam",
                out_dir=config['out_dir'],
                sample_id=config["reads"].keys(),
                param_group=config['bowtie2_params_dict'].keys())


# rule all: #Initial Rule in snakemake
#     input:
#         "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns",
#         "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns",
#         "{0}".format(config['out_dir'])+"/last_trained_alignments/last_counts_alns"
#
#         #expand("{out_dir}/bowtie2_alignments/{sample_id}.{params}/alns.sam",
#         #         out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=bowtie2_paramspace.instance_patterns)
#         # expand("{out_dir}/last_alignments/{sample_id}/{params}/alns.tab",
#         #          out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=last_paramspace.instance_patterns)
#         #expand("{out_dir}/last_alignments/{sample_id}.tab",
#         #         out_dir = config["out_dir"],sample_id=config["reads"].keys())
# rule count_last_trained_alns:
#     input:
#         aln_files = expand("{out_dir}/last_trained_alignments/{sample_id}.tab",
#                            out_dir = config["out_dir"],sample_id=config["reads"].keys())
#     output:
#         "{0}".format(config['out_dir'])+"/last_trained_alignments/last_counts_alns"
#     run:
#         out = open(output[0],"w")
#         out.write("samples\t"+"\taln_counts\n")
#         for aln_filename in input.aln_files:
#             with open(aln_filename) as f:
#                 aln_count = 0
#                 for line in f:
#                     if not line.startswith("#"):
#                         aln_count +=1
#
#             #separate aln_filename into the wildcards and params
#             filename_split = aln_filename.split("/")
#             index = filename_split.index("last_trained_alignments")
#             sample_id = filename_split[index+1]
#             towrite = sample_id+"\t"+f"{aln_count}"+"\n"
#             out.write(towrite)
#         out.close()
#
# rule count_last_alns:
#     input:
#         aln_files=expand("{out_dir}/last_alignments/{sample_id}/{params}/alns.tab",
#                          out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=last_paramspace.instance_patterns)
#     output:
#         "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns"
#     run:
#         out = open(output[0],"w")
#         param_names = list(last_params_pd)
#         out.write("samples\t"+".".join(param_names)+"\taln_counts\n")
#         for aln_filename in input.aln_files:
#             with open(aln_filename) as f:
#                 aln_count = 0
#                 for line in f:
#                     if not line.startswith("#"):
#                         aln_count +=1
#
#             #separate aln_filename into the wildcards and params
#             filename_split = aln_filename.split("/")
#             index = filename_split.index("last_alignments")
#             sample_id = filename_split[index+1]
#             params_value=filename_split[index+2:-1]
#             values = ".".join([ pair.split("~")[-1] for pair in params_value])
#             towrite = sample_id+"\t"+values+"\t"+f"{aln_count}"+"\n"
#             out.write(towrite)
#         out.close()
#
#
#
# rule count_bowtie2_alns:
#     input:
#         aln_files= expand("{out_dir}/bowtie2_alignments/{sample_id}/{params}/alns.sam",
#                            out_dir = config["out_dir"],sample_id=config["reads"].keys(),params=bowtie2_paramspace.instance_patterns)
#     output:
#         "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns"
#     run:
#         out = open(output[0],"w")
#         param_names = list(bowtie2_params_pd)
#         out.write("samples\t"+".".join(param_names)+"\taln_counts\n")
#         for aln_filename in input.aln_files:
#             with open(aln_filename) as f:
#                 aln_count = 0
#                 for line in f:
#                     if not line.startswith("#"):
#                         aln_count +=1
#
#             #separate aln_filename into the wildcards and params
#             filename_split = aln_filename.split("/")
#             index = filename_split.index("bowtie2_alignments")
#             sample_id = filename_split[index+1]
#             params_value=filename_split[index+2:-1]
#             values = ".".join([ pair.split("~")[-1] for pair in params_value])
#             towrite = sample_id+"\t"+values+"\t"+f"{aln_count}"+"\n"
#             out.write(towrite)
#         out.close()
