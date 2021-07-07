# Introduction
This is a pipeline for drug-resistance profiling from HIV whole genome NGS data.
It takes in a HIV reference genome (currently HXB2) and a dataset of reads in fastq format, 
and does the following: aligns the reads to the reference, performs variant calling, and queries the HIVDB sysystem for the presence and degree of drug resistance.


# Installation

Prepare a conda environment using the environment.yaml file inside the config directory.
TODO: incorporate conda inside snakemake so that the only requirements are conda and snakemake.

# Running

With the conda environment activated run:
```
snakemake --configfile config/config.drm.yaml -np 
```
to do a dry-run. If everything seems ok, then run:
```
snakemake --configfile config/config.drm.yaml --cores all
```
