# R_utilities #

Last version: 17-07-2020

This repository contains scripts to make repetitive task that can be recycled across multiple projects.
Although many scripts are related to genomics (the field in which I work) others may be used for non-genomic-related purposes. 

I will be adding the scripts once they are documented and tested.


## Installation

Make sure you have installed:

* R (version >= 3.6.1)
* All scripts use *optparse* R package, be sure you have it installed.


### Extract genes

1. Get all the protein coding genes annotations in human genome from ENSEMBL

Requires:

  - biomaRt (Biocondcutor package)

Given a human genome version [hg19|hg18] returns the list of all protein coding genes, with their IDs and genomic coordinates.

Parameters:

    *-g* : human genome version [hg19|hg18]
    *-o* : output folder

Example:
```
Rscript R-scripts/List_of_ENSEMBL_protein_coding_genes.R -o Results_dir -g hg19
```