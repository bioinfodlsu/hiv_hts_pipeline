# rule mark_duplicates:
#     input:
#         "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam",
#     output:
#         "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.marked.bam",
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         "gatk MarkDuplicatesSpark -I {input} -O {output}"
#
# rule recal_gatk:
#     input:
#
#     output:
#
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         "gatk MarkDuplicatesSpark -I {input} -O {output}"
#
# rule indels_lofreq:
#     input:
#         reference = config["reference"]
#     output:
#         bam = "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.sorted.bam",
#         recal = "{0}".format(config['out_dir'])+"/bowtie2_alignments/{sample_id}/paramgroup_{param_group}/alns.baserecal.bam"
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/lofreq.yaml"
#     shell:
#         "lofreq indelqual --dindel -f {input.reference} -o {output.recal} {output.bam}"
#         """
#         lofreq indelqual -f {input.reference} -b -
#         lofreq alnqual -f $reffa  -b - -o
#         """
