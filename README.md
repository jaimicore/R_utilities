# R_utilities #

Last version: 21-07-2020

This repository contains scripts to make repetitive task that can be recycled across multiple projects.
Although many scripts are related to genomics (the field in which I work) others may be used for non-genomic-related purposes. 

I will be adding the scripts once they are documented and tested.


## Installation

Make sure you have installed:

* R (version >= 3.6.1)
* All scripts use *optparse* R package, be sure you have it installed.


_____
_____


### Extract human protein-coding genes annotations from ENSEMBL

Requires:

  - biomaRt (Biocondcutor)

Given a human genome version [hg19|hg18] returns the list of all protein coding genes, with their IDs and genomic coordinates.

Parameters:

    -g : human genome version [hg19|hg18]
    -o : output folder

Example:
```
Rscript R-scripts/List_of_ENSEMBL_protein_coding_genes.R -o Results_dir -g hg19
```


_____
_____


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


_____
_____


### Generate random coordinates from a VCF or BED file

Requires:

  - data.table (CRAN)
  - GenomicRanges (Bioconductor)

Script to generate random mutations from a BED/VCF template file. 
The random mutations are generated from a uniform distribution, determined from a min and max distance relative to the original position.

In case of required, this script can be used to generate new IDs using a prefix.

Example of input files:

```unix
## Example input VCF file
CHROMOSOME  POSITION  REF ALT SAMPLE
chr1	1671	C	T	Sample1
chr1	5372	G	T	Sample2
chr1	5485	C	T	Sample1
chr1	7729	T	A	Sample3
chr1	11924	G	A	Sample1
```


```unix
## Example input BED file
CHROMOSOME  START END ELEMENT
chr1    10144   10257   Enhancer
chr1    10176   10192   TFBS
chr1    10566   10677   Enhancer
chr1    15217   15230   TFBS
chr1    15627   15645   TFBS
```

Parameters:

    -o : (--output_directory) Output folder
    -t : (--template) Input BED/VCF file
    -f : (--format) Input file extension [BED|VCF]
    -s : (--side) The side where the new random positions will be generated [left|right|both] relative to the original position
    -m : (--max_distance) The maximum number of nucleotides from the original positions where the new random positions will be generated [Default: 10]
    -n : (--min_distance) The minimum number of nucleotides from the original positions where the new random positions will be generated [Default: 1]
    -i : (--new_ids) Column number containing the IDs. When this option is > 0, it will generate new IDs with random numbers. [Default: 0]
    -p : (--id_prefix) Prefix for the new IDs. the new IDs will contain the prefix followed by a random number. E.g. Sample001 [Default: Sample]


Example:
```
## BED file
Rscript R-scripts/Generate_random_coordinates.R -o BED_example -t Coordinates.bed -f bed -s both -m 15 -n 2

## VCF file
Rscript R-scripts/Generate_random_coordinates.R -o VCF_example -t Coordinates.vcf -f vcf -s both -m 10 -n 2 -i 5 -p DONOR
```

**BED example**: the new coordinates will be generated within a window of 2 to 15 nucleotides from the original coordinates.

**VCF example**: the new coordinates will be generated within a window of 2 to 10 nucleotides from the original coordinates.
Note that the new coordinates will not overlap the original ones. 
In addition, the 5th column contains the IDs, and will be replaced using the prefix 'DONOR', e.g., DONOR001.


_____
_____


