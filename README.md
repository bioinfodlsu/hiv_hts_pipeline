# Introduction
This is a pipeline for drug-resistance profiling from HIV whole genome NGS data.
It takes in a HIV reference genome (currently HXB2) and a dataset of reads in (gzipped) fastq format, 
and performs the following steps: 

1. aligns the reads to the reference, 
2. performs variant calling, and
3. queries the HIVDB system for the presence and degree of drug resistance.

Currently, the pipeline uses Bowtie2 for Step1, Lofreq for Step 2, and sierrapy for Step 3.
It has only been tested upto Step 2.

# Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
All other dependencies are handled automatically by Snakemake.

### Install Conda 
Download Miniconda3  installer for Linux from  [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
Installation instructions are [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Once installation is complete, you can test your Miniconda installation by running:
```
$ conda list
```

### Install Snakemake
Snakemake recommends installation via Conda:
```
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```

### Download the pipeline
Clone or download SAMAR from the online  [repository](https://bitbucket.org/hiv_hts_pipeline/hiv_pipeline/).

# Quickstart Guide
Let's try running the pipeline on sample data provided in the `test_data` folder.
With the snakemake conda environment activated, and from the top-level directory (i.e. the one that contains this readme file), run:
```
snakemake --use-conda --configfile config/config.sample.yaml -np
```
to do a dry-run. If snakemake does not complain and everything seems ok, then run:
```
snakemake --use-conda --configfile config/config.sample.yaml --cores all
```
The results can be found inside the newly created directory called `test_result`.

# Running the pipeline on your own data
To the run the pipeline on your own data, you need to specify in a config file the paths to the input data (reads and reference), path to the output directory. Optionally, in this config file, you can also set parameters for the various tools that make up this pipeline. You can use `config/config.drm.yaml` as a template.  Once the configfile is ready, run the pipeline like above:
```
snakemake --use-conda --configfile config/config.sample.yaml -np
```
for a dry run, and 
```
snakemake --use-conda --configfile config/config.sample.yaml --cores all
```
for the actual run.
