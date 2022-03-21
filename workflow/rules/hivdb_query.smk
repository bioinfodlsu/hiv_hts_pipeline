rule prepare_query:
    input:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants_{aligner}.vcf"
    output:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants_{aligner}.json"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/vcfToHivdb.py --vcfFile {input} --outFile {output}"

rule tabulate_results:
    input:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants_{aligner}.json"
    output:
        "{0}".format(config['out_dir'])+"/drug_resistance_report/{sample_id}/paramgroup_{param_group}/results_{aligner}.txt"
    conda:
        "../envs/vcfToHivdb.yaml"
    shell:
        "python workflow/scripts/tabulateJson.py --jsonFile {input} --outFile {output}"
