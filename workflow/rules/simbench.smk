rule simulate_all:
    input:
        expand(
            f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bai",
            mutrate=[f"mr_{x:.3f}" for x in [i/1000 for i in range(1, 21)]]
        )
    
rule simulate_haplotypes:
    input:
        reference="test_data/refseqs_cleaned/K03455.fasta"
    params:
        simbench="python /Users/gcoe/Documents/GitHub/vp-analysis/V-pipe/workflow/scripts/simBench.py",
        vpipe_dir="/Users/gcoe/Documents/GitHub/vp-analysis/V-pipe/my_haplotypes",
        outdir=lambda wildcards: f"test_data/simulated_reads/my_haplotypes/{wildcards.mutrate}",
        num_haplotypes=5,
        mutation_rate=lambda wildcards: float(wildcards.mutrate.replace("mr_", "")),
        seed=42,
    output:
        touch("test_data/simulated_reads/my_haplotypes/{mutrate}/haplotypes.done")
    conda:
        "../envs/simbench.yaml"
    shell:
        """
        {params.simbench} -f {input.reference} \
            -n {params.num_haplotypes} --tree-like \
            -mr {params.mutation_rate} -s {params.seed} \
            -oh {params.vpipe_dir} -o haplotypes \
            > >(tee {params.vpipe_dir}/simbench.out.log) 2>&1

        mkdir -p {params.outdir}
        mv {params.vpipe_dir}/*.fasta {params.outdir}
        mv {params.vpipe_dir}/*.log {params.outdir}
        """

rule simulate_reads:
    input:
        done="test_data/simulated_reads/my_haplotypes/{mutrate}/haplotypes.done"
    params:
        simbench="python /Users/gcoe/Documents/GitHub/vp-analysis/V-pipe/workflow/scripts/simBench.py",
        outdir=lambda wildcards: f"test_data/simulated_reads/my_reads/{wildcards.mutrate}",
        haplotypes_dir=lambda wildcards: f"test_data/simulated_reads/my_haplotypes/{wildcards.mutrate}",
        num_haplotypes=5,
        seed=42,
    output:
        touch("test_data/simulated_reads/my_reads/{mutrate}/simreads_R.done")
    conda:
        "../envs/simbench.yaml"
    shell:
        """
        {params.simbench} -n {params.num_haplotypes} -q -art art_illumina -v \
            -oh {params.haplotypes_dir} -or {params.outdir} -o reads \
            > >(tee {params.outdir}/simbench.out.log) 2>&1

        mkdir -p {params.outdir}/reads
        mv {params.outdir}/*.fastq {params.outdir}/reads || true
        mv {params.outdir}/*.fq-e {params.outdir} || true
        mv {params.outdir}/*.sam {params.outdir} || true
        mv {params.outdir}/*.aln {params.outdir} || true
        mv {params.outdir}/*.log {params.outdir} || true
        """

rule recalibrate_simulated_sam:
    input:
        done="test_data/simulated_reads/my_reads/{mutrate}/simreads_R.done",
    params:
        script="workflow/scripts/recalibrate_simulated_sam.py",
        fasta=lambda wildcards: f"test_data/simulated_reads/my_haplotypes/{wildcards.mutrate}/haplotypes.fasta",
        sam_file=lambda wildcards: f"test_data/simulated_reads/my_reads/{wildcards.mutrate}/simreads_R.sam",
        out_sam=lambda wildcards: f"test_data/simulated_reads/my_reads/{wildcards.mutrate}/simreads_R.recalibrated.sam",
    output:
        done="test_data/simulated_reads/my_reads/{mutrate}/simreads_R.recalibrated.done",
    conda:
        "../envs/simbench.yaml"
    shell:
        """
        python {params.script} --haplotypes {params.fasta} --input {params.sam_file} --output {params.out_sam}
        touch {output.done}
        """

rule sam_to_bam:
    input:
        done="test_data/simulated_reads/my_reads/{mutrate}/simreads_R.recalibrated.done",
    params:
        outdir=lambda wildcards: f"test_data/simulated_reads/my_reads/{wildcards.mutrate}",
        sam_file=lambda wildcards: f"test_data/simulated_reads/my_reads/{wildcards.mutrate}/simreads_R.recalibrated.sam",
    output:
        bam=f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bam",
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        samtools view -S -bh {params.sam_file} | samtools sort -o {output.bam}
        """

rule index_bam:
    input:
        bam=f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bam",
    output:
        bai=f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bai",
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule lofreq_simulated:
    input:
        reference = f"test_data/references/hxb2.fasta",
        bam = f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bam",
        bai = f"test_data/simulated_reads/my_reads/{{mutrate}}/simreads_R.recalibrated.sorted.bai",
    output:
        f"test_data/simulated_reads/my_reads/{{mutrate}}/variants/simreads_R.vcf",
    params:
        outdir=f"test_data/simulated_reads/my_reads/{{mutrate}}/variants",
    log:
        f"test_data/simulated_reads/my_reads/{{mutrate}}/variants/simreads_R.lofreq.log",
    conda:
        "../envs/linear.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir} &&
        lofreq faidx {input.reference} &&
        lofreq call --no-default-filter -N -f {input.reference} -o {output} {input.bam} 2> {log}
        """

rule vcf_to_aavf_simulated:
    input:
        f"test_data/simulated_reads/my_reads/{{mutrate}}/variants/simreads_R.vcf"
    output:
        f"test_data/simulated_reads/my_reads/{{mutrate}}/variants/simreads_R.aavf" 
    params:
        out_dir = f"test_data/simulated_reads/my_reads/{{mutrate}}/variants"
    conda:
        "../envs/vcf_to_hivdb.yaml"
    log:
        f"test_data/simulated_reads/my_reads/{{mutrate}}/variants/vcf_to_aavf_simreads_R.log"
    shell:
        """
        python workflow/scripts/vcf_to_aavf_2.py --vcf {input} --aavf {output} 2> {log}
        """
