rule prepare_query:
    input:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants.vcf"

    output:
        "{0}".format(config['out_dir'])+"/variants/{sample_id}/paramgroup_{param_group}/variants.sierra"

    script:
        vcfToHivdb.py {input} {output}
