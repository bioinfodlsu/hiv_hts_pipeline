import os


def get_read1(wildcards):
    return "{}/filtered_reads/{}/{}_1_subsampled.fq".format(config['out_dir'],wildcards.sample_id,wildcards.sample_id);

def get_read2(wildcards):
    return "{}/filtered_reads/{}/{}_2_subsampled.fq".format(config['out_dir'],wildcards.sample_id,wildcards.sample_id);
    
def get_reads_count(wildcards):
    read_counts = {}
    read1 = get_read1(wildcards)
    read2 = get_read2(wildcards)

    for a in reads1:
        read_counts[a] = get_read_count(read1,read2)

    a = os.popen("echo $(zcat " + read1 +"|wc -l)/4|bc")
    a = int(a.read())
    b = os.popen("echo $(zcat " + read2 +"|wc -l)/4|bc")
    b = int(b.read())

    return a+b

rule count_all:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/to_{reference_name}/last_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/last_trained_alignments/to_{reference_name}/last_trained_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/to_{reference_name}/bowtie2_counts_alns.tsv"

rule count_last_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               reference_name = config['reference_name'],
               param_group = config['aligner_params_dict'].keys(),
               country = config['country']
              )
    output:
        "{0}".format(config['out_dir'])+"/last_alignments/to_{reference_name}/last_counts_alns.tsv"
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule count_last_trained_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               reference_name = config['reference_name'],
               param_group = config['aligner_params_dict'].keys(),
               country = config['country']
              )
    output:
        "{0}".format(config['out_dir'])+"/last_trained_alignments/to_{reference_name}/last_trained_counts_alns.tsv"
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule count_bowtie2_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               reference_name = config['reference_name'],
               param_group = config['aligner_params_dict'].keys(),
               country = config['country']
              )
    output:
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/to_{reference_name}/bowtie2_counts_alns.tsv"
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule plot_stats:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/to_{reference_name}/last_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/last_trained_alignments/to_{reference_name}/last_trained_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/to_{reference_name}/bowtie2_counts_alns.tsv"
    output:
        "{0}".format(config['out_dir'])+"/alignment_analysis/to_{reference_name}/alignments_plot_{country}.pdf"
    params:
        counts = "{0}".format(config['out_dir'])+"/alignment_analysis/to_{reference_name}/counts_{country}.tsv"
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/plot_stats.py --tsvFile {input} --outFile {output} --countFile {params.counts}"

# rule compare_alns:
#     input:
#         expand("{out_dir}/bowtie2_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
#                out_dir = config['out_dir'],
#                sample_id =  config['reads'].keys(),
#                reference_name = config['reference_name'],
#                param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
#                country = config['country']
#               ),
#         expand("{out_dir}/last_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
#                out_dir = config['out_dir'],
#                sample_id =  config['reads'].keys(),
#                reference_name = config['reference_name'],
#                param_group = config['last_params_dict'].keys(),
#                country = config['country']
#               ),
#         expand("{out_dir}/last_trained_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
#                out_dir = config['out_dir'],
#                sample_id =  config['reads'].keys(),
#                reference_name = config['reference_name'],
#                param_group = config['last_params_dict'].keys(),
#                country = config['country']
#               ),
#         expand("{out_dir}/blast_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/alns.sam",
#                out_dir = config['out_dir'],
#                sample_id =  config['reads'].keys(),
#                reference_name = config['reference_name'],
#                param_group = config['last_params_dict'].keys(),,
#                country = config['country']
#               )
#     output:
#         "{0}".format(config['out_dir'])+"/alignment_analysis/to_{reference_name}/plot_alns_{country}.tsv"
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/blast.yaml"
#     shell:
#         "python workflow/scripts/blast_alns.py --inFile {input} --outFile {output}"

# rule plot_alns:
#     input:
#         "{0}".format(config['out_dir'])+"/alignment_analysis/to_{reference_name}/plot_alns.tsv"
#     output:
#         "{0}".format(config['out_dir'])+"/alignment_analysis/to_{reference_name}/alignments_comparison.png"
#     conda:
#         "../envs/count.yaml"
#     shell:
#         "python workflow/scripts/plot_alns.py --tsvFile {input} --outFile {output}"
