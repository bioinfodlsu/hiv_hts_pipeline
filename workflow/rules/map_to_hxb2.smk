rule p2_last_index:
    input:
        config["hxb2"]
    output:
        touch("{0}/phase_two_index/index.done".format(config["out_dir"]))
    params:
        index_basename = "{0}".format(config['out_dir'])+"/phase_two_index/hxb2_index"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        lastdb -P0 -uNEAR -R01 {params.index_basename} {input}
        """

rule p2_last_train:
    input:
        subtype = config["reference"],
        flag = "{0}/phase_two_index/index.done".format(config["out_dir"])
    output:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/hxb2_{reference_name}.mat"
    params:
        index_basename = "{0}".format(config['out_dir'])+"/phase_two_index/hxb2_index"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        last-train -P0 --verbose --revsym --matsym --gapsym -E0.05 -C2 {params.index_basename} {input.subtype} > {output}
        """

rule map_manyto1:
    input:
        subtype = config["reference"],
        freq = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/hxb2_{reference_name}.mat"
    output:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping_1.maf"
    params:
        index_basename = "{0}".format(config['out_dir'])+"/phase_two_index/hxb2_index"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        lastal -m50 -E0.05 -C2 -p {input.freq} {params.index_basename} {input.subtype} | 
        last-split -m1 > {output}
        """

rule map_1to1:
    input:
        manyto1 = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping_1.maf"
    output:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.maf"
    params:
        subtype = config['reference_name']
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        maf-swap {input.manyto1} |
        awk '/^s/ {{$2 = (++s % 2 ? "{params.subtype}." : "hxb2.") $2}} 1' |
        last-split -m1 |
        maf-swap > {output}
        """

rule maf_to_sam:
    input:
        query = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.maf",
        dict = "{0}".format(config["out_dir"])+"/filtered_reads/{sample_id}/{sample_id}.dict"
    output:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.sam"
    params:
        subtype = config['reference_name']
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        "maf-convert sam {input.query} -f {input.dict} > {output}"

rule filter_alignments:
    input:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.maf"
    output:
        "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.tab"
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/last.yaml"
    shell:
        """
        last-postmask {input} |
        maf-convert -n tab |
        awk -F'=' '$2 <= 1e-5' > {output}
        """

# rule p2_prepare_query:
#     input:
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf",
#         maf = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.maf"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.json"
#     conda:
#         "../envs/vcfToHivdb.yaml"
#     shell:
#         "python workflow/scripts/ref_to_hivdb.py --vcfFile {input.vcf} --mafFile {input.maf} --outFile {output}"


# rule p2_prepare_query_e:
#     input:
#         vcf = "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}.vcf",
#         maf = "{0}".format(config['out_dir'])+"/phase_two_subtype_alignments/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/mapping.maf"
#     output:
#         "{0}".format(config['out_dir'])+"/variants/{country}/{sample_id}_to_{reference_name}/paramgroup_{param_group}/variants_{aligner}_{lofreq_param_group}.json"
#     conda:
#         "../envs/vcfToHivdb.yaml"
#     shell:
#         "python workflow/scripts/ref_to_hivdb.py --vcfFile {input.vcf} --qual  --mafFile {input.maf} --outFile {output}"

# rule maptohxb2:
#     input:
#         hxb2 = config["hxb2"]
#         subtype = config["reference"]
#     output:
#         "{0}".format(config['out_dir'])+"/phase_two/{country}/{sample_id}/paramgroup_{param_group}/index"
#     params:
    
#     threads:
#         workflow.cores/len(config["reads"])
#     conda:
#         "../envs/last.yaml"
#     shell:
#         """
#         lastdb -P8 -uMAM8 myDB genome1.fa
#         last-train -P8 --revsym -D1e9 --sample-number=5000 myDB genome2.fa > my.train
#         lastal -P8 -D1e9 -m100 --split-f=MAF+ -p my.train myDB genome2.fa > many-to-one.maf
#         last-split -r many-to-one.maf | last-postmask > out.maf
#         """
