# rule prepare_query_qual:
#     input:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.json"
#     params:
#         QUAL = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]['QUAL']
#     conda:
#         "../envs/vcfToHivdb.yaml"
#     shell:
#         "python workflow/scripts/vcfToHivdb.py --vcfFile {input} --qual {params.QUAL} --outFile {output}"

rule prepare_query:
    input:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.json"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/vcfToHivdb.py --vcfFile {input} --outFile {output}"


rule tabulate_results:
    input:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.json"
    output:
        "{0}".format(config['out_dir'])+"/drug_resistance_report/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/results_{aligner}.txt"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/tabulateJson.py --jsonFile {input} --outFile {output}"

# rule tabulate_results_e:
#     input:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.json"
#     output:
#         "{0}".format(config['out_dir'])+"/drug_resistance_report/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/results_{aligner}_{lofreq_param_group}.txt"
#     conda:
#         "../envs/vcfToHivdb.yaml"
#     shell:
#         "python workflow/scripts/tabulateJson.py --jsonFile {input} --outFile {output}"
