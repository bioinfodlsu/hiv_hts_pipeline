rule count_all:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/last_trained_alignments/last_trained_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns.tsv"

rule count_last_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              )
    output:
        temp("{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns.tsv")
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule count_last_trained_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              )
    output:
        temp("{0}".format(config['out_dir'])+"/last_trained_alignments/last_trained_counts_alns.tsv")
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule count_bowtie2_alns:
    input:
        expand("{0}".format(config['out_dir'])+"/bowtie2_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              )
    output:
        temp("{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns.tsv")
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/count_alns.py --inFile {input} --outFile {output}"

rule plot_stats:
    input:
        "{0}".format(config['out_dir'])+"/last_alignments/last_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/last_trained_alignments/last_trained_counts_alns.tsv",
        "{0}".format(config['out_dir'])+"/bowtie2_alignments/bowtie2_counts_alns.tsv"
    output:
        "{0}".format(config['out_dir'])+"/alignment_analysis/alignments_plot_{country}.pdf"
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/plot_stats.py --tsvFile {input} --outFile {output}"

rule compare_alns:
    input:
        expand("{out_dir}/bowtie2_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              ),
        expand("{out_dir}/last_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              ),
        expand("{out_dir}/last_trained_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              ),
        expand("{out_dir}/blast_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam",
               out_dir = config['out_dir'],
               sample_id =  config['reads'].keys(),
               param_group = config['{}_params_dict'.format(config['aligner'])].keys(),
               country = config['country']
              )
    output:
        "{out_dir}/alignment_analysis/plot_alns.tsv"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/blast.yaml"
    shell:
        "python workflow/scripts/blast_alns.py --inFile {input} --outFile {output}"

rule plot_alns:
    input:
        "{out_dir}/alignment_analysis/plot_alns.tsv"
    output:
        "{out_dir}/alignment_analysis/alignments_comparison.png"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/count.yaml"
    shell:
        "python workflow/scripts/plot_alns.py --tsvFile {input} --outFile {output}"
