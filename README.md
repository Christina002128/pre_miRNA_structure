# Tropic Associate Bioinformatician Assessment - miRNA structure predictions


## Introduction

Welcome to your assignment for the role of Associate Bioinformatician. This assignment will be available to you until 4th December 2023 at 5pm and is composed of 2 parts:

- Four theoretical questions, detailed below. Please choose **two** to answer as separate issues.
- A coding exercise to complete a Snakemake pipeline to extract miRNA precursor sequences from a genome and predict their secondary structure.


This assignment is designed to assess your technical skills as well as your problem solving and communication abilities. Therefore, communication with the team around your progress is crucial, **there are no stupid questions**.

## Theoretical questions

For each question, please aim to give high-level answers in a few sentences only. You are free to use assumptions in your answers, however please make sure they are clearly stated in your answers.

1. One of Tropic's main activity revolves sequencing small RNA for analysis of miRNA expression in different tissues. Please outline what key challenges can be encountered
when handling sRNAseq data.

2. Following identification of a novel non-transgenic, gene edited banana plant, we aim to submit this candidate to a regulatory body for exemption from regulation. To this end, we are required to demonstrate the complete absence of transgenic sequences derived from our vector. To this extent, whole genome sequencing was performed on the plant of interest. Please describe the process you would follow from obtaining the raw sequencing data to proving absence of transgenic material in the genome.

3. Whilst scouring the literature, you found a gene of interest providing drought resistance in Coconut (_Cocos nucifera_). You want to know if a similar gene is present in the coffee genome (_Coffea arabica_). Please detail the steps you would undertake.

4. GEiGS Biocompute (https://www.geigs.com/) predicts solutions allowing to redirect endogenous miRNAs against novel endogenous or pathogenic targets. Please explain potential challenges you might encounter during implementation of GEiGS solutions _in planta_.


## miRNA structure predictions exercise

You are in charge of leading a project to accurately predict miRNA structures for the GEiGS Biocompute platform. Unfortunately, one of your colleagues just won the lottery and abruptly left the team, leaving the pipeline unfinished. Your first task is to complete this Snakemake pipeline to extract miRNA precursor sequences from a genome of your choice and predict their secondary structure. You're not alone in this project and can count on the help and collaboration of your team.

## Instructions

- You have "developer level" access to this repository, which means you won't be able
  to push changes to the main branch. You'll have to branch out from main and
  submit merge requests for review.
- Feel free to do a single merge request for the whole thing, or several where
each makes an incremental change
- Look for hints in source code comments and complete the pipeline.
- Use the merge request's comments section (and any other available features
available in GitLab that you may find useful) to keep the team updated on your progress.

### Input files

You'll need to download a reference genome and its corresponding GFF3 file from
a public database (e.g. ENSEMBL), and update the `config/config.yaml` file
accordingly.

### Output files

The final output files should contain 4 tab-separated columns:
- The miRNA ID (any identifier is fine)
- The miRNA precursor sequence
- The miRNA precursor structure in dot-bracket format
- The structure's MFE

### Hints

#### Software requirements

You will need conda to manage the required software packages for the pipeline. We suggest you install conda via
[mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

After this, you should set up bioconda to access a vast library of
bioinformatics software by running these commands:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Once you have conda installed and activated you can install Snakemake with
`mamba install snakemake`.

#### Running the pipeline

As it is, by completing the TODO item in line 14 you will have a running
pipeline which you can execute with `snakemake --use-conda -c`. This will create
dummy output. We suggest you make sure you can obtain all dummy files before you start coding
your solutions.

#### Software dependencies

You should add any required software to each rule's environment yaml file,
located in `workflow/envs`. As an example, the environments currently
install the latest version of the `bedtools` package
(which is not necessarily needed depending on how you go about implementing each
rule).

#### Input files 

- Make sure you test the pipeline with at least 2 organisms (at the same time)
- It's a good idea to choose smaller organisms to keep things quick

#### Processing GFF files

You may choose writing your own code to parse the GFF files, or resort to a
library (e.g. bcbio).

#### Predicting secondary structure

There are several libraries that can predict secondary RNA structure from sequence, one example is ViennaRNA.
