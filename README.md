# Pangenome Graph-Based HIV Drug Resistance Profiling Pipeline

A Snakemake-based workflow for detecting HIV drug resistance mutations (DRMs) using a pangenome graph reference built from diverse HIV-1 group M genomes. This pipeline leverages the **VG toolkit** and related tools to construct variation graphs, align NGS reads, call variants, and convert results into AAVF format for drug resistance mutation identification and interpretation.


## What This Pipeline Does

This pipeline performs the following key steps:

1. **Input Preparation**

   * Accepts a curated set of HIV reference genomes (FASTA) and NGS reads (FASTQ).

2. **Graph Construction**

   * Performs multiple sequence alignment using `MAFFT`.
   * Constructs a variation graph using `vg construct`.

3. **Graph Indexing and Read Alignment**

   * Indexes the graph using `vg index`, `vg gbwt`, and related tools.
   * Aligns reads to the graph using `vg giraffe` or `vg map`.

4. **Variant Calling**

   * Calls variants from graph alignments using `vg call`.
   * Generates a VCF file per sample.

5. **Drug Resistance Mutation Identification and Interpretation**

   * Converts VCF to AAVF using a custom Python script.
   * Queries the AAVF file against **HIVdb** using `sierrapy`.
   * Outputs an HIV drug resistance profile.

7. **Benchmarking and Evaluation**

   * Compares alignment rate, accuracy, variant detection, runtime, and memory usage against a traditional linear-reference pipeline.


## Installation

This pipeline requires:

* [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* [Snakemake](https://snakemake.readthedocs.io)

All tool-specific dependencies are handled via Snakemake’s `--use-conda` feature.

### Step 1: Install Miniconda

Follow the [official instructions here](https://docs.conda.io/en/latest/miniconda.html) to install Miniconda3 for your platform.

Verify the installation:

```bash
conda list
```

### Step 2: Install Snakemake

We recommend using **mamba** for faster environment setup:

```bash
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n \<environment_name\> snakemake
conda activate \<environment_name\>
```

Test the installation:

```bash
snakemake --help
```


## Download the Pipeline

Clone this repository:

```bash
git clone https://github.com/kimileeee/pangenomics-hiv.git
cd pangenomics-hiv
```


## Quick Start Guide

### Dry Run

To preview the workflow and see what will run:

```bash
snakemake --use-conda --cores N -np
```

Replace `N` with the number of cores you wish to use.

### Run the Pipeline

To execute the pipeline:

```bash
snakemake --use-conda --cores N -p
```

### Visualize DAG (Directed Acyclic Graph)

To generate a PDF of the pipeline's rule structure:

```bash
snakemake --use-conda --cores all --dag | dot -Tpdf > dag.pdf
```


## Running the Pipeline on Your Own Data

### 1. Configure `config/config.yaml`

Edit the `config.yaml` file to set:

* Paths to your **reference FASTA** and **FASTQ** files
* Output directory
* Optionally set `linear: true` to run the linear pipeline instead of the graph-based one
* Optionally set `use_clustered_refs: true` to enable reference clustering before graph construction
* Graph alignment method (`giraffe` or `map`)

Make sure all paths are correct and that FASTA/FASTQ files are in standard format. You can leave unused samples or references commented out.


## Citation & References



## Contact

For questions, bug reports, or contributions, feel free to [open an issue](https://github.com/kimileeee/pangenomics-hiv/issues).

