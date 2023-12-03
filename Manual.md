# README.md for Structure Prediction of miRNA Precursor Pipeline

------------------------------------------------------------------------

## Overview

This pipeline is designed for the extraction, analysis, and conversion of miRNA precursor sequences and structures. Leveraging the Snakemake workflow management system, it automates the process of precursors sequence extraction, structure prediction, and data formatting, facilitating efficient and reproducible analyses.

------------------------------------------------------------------------

## Scripts and Environment Files

-   **Snakefile**: Defines the workflow rules and dependencies.
-   **extract_precursors.py**: Extracts miRNA precursor sequences.
-   **convert_to_tsv.py**: Converts data into TSV format.
-   **config/config.yaml**: Conda configuration for Snakefile.
-   **extract_precursors.yaml**: Conda environment for `extract_precursors.py`.
-   **predict_structure.yaml**: Conda environment for RNA structure prediction.

------------------------------------------------------------------------

## Dependencies

-   [Snakemake](https://snakemake.readthedocs.io)
-   [BioPython](https://biopython.org)
-   [BCBio-GFF](https://github.com/chapmanb/bcbb/tree/master/gff)
-   [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)

------------------------------------------------------------------------

## Installation

1.  Install Snakemake following the [official guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

2.  Clone this repository:

    ``` bash
    git clone https://gitlab.com/tropic-assignment-cz/mirna-str-predict.git
    cd mirna-str-predict
    ```

3.  Create and activate the conda environments:

    ``` bash
    conda env create -f workflow/envs/extract_precursors.yaml
    conda env create -f workflow/envs/predict_structure.yaml 
    ```

------------------------------------------------------------------------

## Usage

1.  Download the genome sequence file in fasta format.

2.  Download the corresponding annotation file in GFF format, which annotates precursor sequences as *miRNA_primary_transcript*.

3.  Modify `config/config.yaml` to specify your genome sequence and annotation files, and the number of sequences you want to process.

4.  Run the Snakemake pipeline:

    ``` bash
    snakemake -c <number of cores>
    ```

------------------------------------------------------------------------

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

------------------------------------------------------------------------

## Contact

For questions or support, please contact [[zxy002128\@163.com](mailto:zxy002128@163.com){.email}].

------------------------------------------------------------------------
