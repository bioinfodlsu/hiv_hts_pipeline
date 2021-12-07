rule prepare_query:
    input:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants.vcf"

    output:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants.json"

    conda:
        "../envs/vcfToHivdb.yaml"

    # script:
    #     "../scripts/vcfToHivdb.py"
    
    shell:
        "python workflow/scripts/vcfToHivdb.py --vcfFile {input} --outFile {output}"


rule tabulate_results:
    input:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants.json"

    output:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/results.txt"

    conda:
        "../envs/vcfToHivdb.yaml"

    # script:
    #     "../scripts/tabulateJson.py"

    shell:
        "python workflow/scripts/tabulateJson.py --jsonFile {input} --outFile {output}"
