rule blast_db:
    input:
        reference = config["reference"]
    output:
        touch("{0}/blast_index/index.done".format(config["out_dir"]))
    params:
        dbtype = "nucl",
        index_basename = "{0}/blast_index/index".format(config["out_dir"])
    threads:
        workflow.cores/len(config["reads"])
    conda:
        "../envs/blast.yaml"
    shell:
        "makeblastdb -in {input.reference} -dbtype {params.dbtype} -out {params.index_basename}"

rule align_blast:
    input:
        reference = config["reference"],
        sub1 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_1_subsampled.fq",
        sub2 = "{0}".format(config['out_dir'])+"/filtered_reads/{sample_id}/{sample_id}_2_subsampled.fq"
    output:
        "{0}".format(config["out_dir"])+"/blast_alignments/{country}/{sample_id}/paramgroup_{param_group}/alns.sam"
    params:
        index_basename = "{0}/blast_index/index".format(config["out_dir"]),
        evalue = 10e-25,
        outfmt = 6
    threads:
        workflow.cores
        # workflow.cores/len(config["reads"])
    conda:
        "../envs/blast.yaml"
    shell:
        """
        fastq-interleave {input.sub1} {input.sub2} |
        seqkit fq2fa |
        blastn -query - -subject {input.reference} -num_threads {threads} -outfmt {params.outfmt} -out {output}
        """