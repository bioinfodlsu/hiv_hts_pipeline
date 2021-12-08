# Introduction
This is a pipeline for drug-resistance profiling from HIV whole genome NGS data.
It takes in a HIV reference genome (currently HXB2) and a dataset of reads in (gzipped) fastq format, 
and performs the following steps: 

1. aligns the reads to the reference, 
2. performs variant calling, and
3. queries the HIVDB system for the presence and degree of drug resistance (requires internet connection).

Currently, the pipeline uses [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for Step1, [Lofreq](https://dx.doi.org/10.1093%2Fnar%2Fgks918) for Step 2, and [sierrapy](https://github.com/hivdb/sierra-client/tree/master/python) for Step 3. 

This code is under development. We hope to test and add more options.

# Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
All other dependencies are handled automatically by Snakemake.

### Install Conda 
Download Miniconda3  installer for Linux from [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers), or for macOS from [here](https://docs.conda.io/en/latest/miniconda.html#macos-installers)
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
Clone this pipeline by clicking the Clone button on the top-right of this page,
or download it by clicking the ellipsis next to the Clone button.

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
Interpretation of drug-resistance as provided by sierrapy can be found inside the `drug_resistance_report` folder.
Intermediate files such as the read-to-reference alignments and variant calls can be found in their respective folders.

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

# Contact
This is an ongoing work. If you have questions, concerns, issues, or suggestions, please contact: 
**Anish Shrestha, 
Bioinformatics Lab, De La Salle University Manila at anish.shrestha@dlsu.edu.ph** .
