# R_utilities #

Last version: 17-07-2020

This repository contains scripts to make repetitive task that can be recycled across multiple projects.
Although many scripts are related to genomics (the field in which I work) others may be used for non-genomic-related purposes. 

I will be adding the scripts once they are documented and tested.


## Installation

Make sure you have installed:

* R (version >= 3.6.1)
* All scripts use *optparse* R package, be sure you have it installed.


### Extract human protein-coding genes annotations from ENSEMBL

Requires:

  - biomaRt (Biocondcutor package)

Given a human genome version [hg19|hg18] returns the list of all protein coding genes, with their IDs and genomic coordinates.

Parameters:

    -g : human genome version [hg19|hg18]
    -o : output folder

Example:
```
Rscript R-scripts/List_of_ENSEMBL_protein_coding_genes.R -o Results_dir -g hg19
```


### Convert tab-delimited file to RData and viceversa

Requires:

  - data.table (CRAN)

Script to convert an RData table into a text file, or a text-file (table) into an RData object.

Parameters:

    -f : (--from) input file extension
    -t : (--to) output file extension
    -i : (--input_file) input file
    -o : (--out_folder) output folder
    -r : (--rownames) Indicates wheter the output text table should have rownames [0|1] [Default: 0]
    -l : (--colnames) Indicates wheter the output text table should have colnames [0|1] [Default: 0]
    -e : (--extension) Extension for the output text table [Default: tab]
    
Example:
```
## From Table to RData
Rscript R-scripts/RData_Tab_converter.R -i TFBSs_hg19.RData -from RData -to tab -o hg19_TFBS -e bed

## From RData to table
Rscript R-scripts/RData_Tab_converter.R -i TFBSs_Unibind_EZ_hg19.bed -from tab -to RData -o hg19_TFBS
```
