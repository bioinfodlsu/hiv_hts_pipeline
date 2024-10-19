rule vcf_to_aavf:
    input:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.aavf"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/vcf_to_aavf.py --vcf {input} --aavf {output}"

rule aavf_filtered:
    input:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.aavf"
    output:
        "{0}".format(config['out_dir'])+"/variants_filtered/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.aavf"
    params:
        filter_param = lambda wildcards: config["lofreq_params_dict"][wildcards.lofreq_param_group]
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/aavf_filter.py --aavf {input} --out {output} -A {params.filter_param[AF]}"
        # "python workflow/scripts/aavf_filter.py --aavf {input} --out {output} -C {params.filter_param[DP]}" 

rule aavf_to_hivdb:
    input:
        "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.aavf"
    output:
        "{0}".format(config['out_dir'])+"/drug_resistance_report/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/results_{aligner}.txt"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/aavf_to_hivdb.py --aavf {input} --out {output}"

rule aavf_to_hivdb_filtered:
    input:
        "{0}".format(config['out_dir'])+"/variants_filtered/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.aavf"
    output:
        "{0}".format(config['out_dir'])+"/drug_resistance_report_filtered/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/results_{aligner}_{lofreq_param_group}.txt"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/aavf_to_hivdb.py --aavf {input} --out {output}"
