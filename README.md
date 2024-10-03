# dataset_collation_tool.py

This repository contains a simple script for assembling datasets used in viral genomic epidemiology. It utilizes the NCBI Datasets command-line interface tool to retrieve and sample relevant genetic diversity.

## Command-Line Arguments
The script accepts the following command-line arguments:

    --taxid                NCBI Taxonomy ID for the target virus.
    --target               Number of sequences to include per location and time period (default = 1).
    --length-filter        Minimum sequence length to be included in the dataset (default = 10,000 nucleotides).
    --ambiguity-filter     Maximum allowed percentage of ambiguous characters (default = 0.25, range: 0 - 1).
    --geo-unit             Geographic level for sampling (default: country) [options: country, continent].
    --time-unit            Temporal level for sampling (default: year) [options: year, year-month].
    --output OUTPUT        Full path for the output directory.
    --tag                  Tag or name for the subsetting scheme.
    --skip-download        Option to skip the download step (assumes the NCBI data package has already been downloaded and all paths are set correctly).

## Example Usage:

        python dataset_collation_tool.py --taxid 64320 --output zika_virus_ncbi_data --ambiguity-filter 0.1 --target 10

In this example, the script will:

* Download all Zika virus data available from NCBI.
* Filter out sequences shorter than 10,000 nucleotides and those with more than 10% ambiguous characters.
* Sample up to 10 sequences per year per country to assemble a globally representative dataset.

The tool currently supports year-month and continent-based sampling. Future versions will offer more flexible sampling schemes.

## Installation and Dependencies

Ensure that all dependencies are installed as specified in the env/dataset_collation_tool.yml file.

## Citation
If you find this script useful, please consider citing this GitHub repository. =)

