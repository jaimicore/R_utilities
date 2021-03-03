# R_utilities #

Last version: 03-03-2021

This repository contains scripts to make repetitive task that can be recycled across multiple projects.
Although many scripts are related to genomics (the field in which I work) others may be used for non-genomic-related purposes. 

I will be adding the scripts once they are documented and tested.

Any question, contact me: j.a.c.mondragon@ncmm.uio.no


## Installation

Make sure you have installed:

* R (version >= 3.6.1)
* All scripts use *optparse* R package, be sure you have it installed.


_____


### Extract human protein-coding genes annotations from ENSEMBL

Requires:

  - biomaRt (Biocondcutor)

Given a human genome version [hg19|hg18] returns the list of all protein coding genes, with their IDs and genomic coordinates.

Parameters:

    -g : Human genome version [hg19|hg18]
    -o : Output folder

Example:
```
Rscript R-scripts/List_of_ENSEMBL_protein_coding_genes.R \
  -o Results_dir                                         \
  -g hg19
```


_____


### Convert tab-delimited file to RData and viceversa

Requires:

  - data.table (CRAN)

Script to convert an RData table into a text file, or a text-file (table) into an RData object.

Parameters:

    -f : (--from)           Input file extension
    -t : (--to)             Output file extension
    -i : (--input_file)     Input file
    -o : (--out_folder)     Output folder
    -r : (--rownames)       Indicates wheter the output text table should have rownames [0|1] [Default: 0]
    -l : (--colnames)       Indicates wheter the output text table should have colnames [0|1] [Default: 0]
    -e : (--extension)      Extension for the output text table [Default: tab]
    
Example:
```
## From Table to RData
Rscript R-scripts/RData_Tab_converter.R \
  -i TFBSs_hg19.RData                   \
  -from RData                           \
  -to tab                               \
  -o hg19_TFBS                          \
  -e bed

## From RData to table
Rscript R-scripts/RData_Tab_converter.R \
  -i TFBSs_Unibind_EZ_hg19.bed          \
  -from tab                             \
  -to RData                             \
  -o hg19_TFBS
```


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

    -o : (--output_directory)   Output folder
    -t : (--template)           Input BED/VCF file
    -f : (--format)             Input file extension [BED|VCF]
    -s : (--side)               The side where the new random positions will be generated [left|right|both] relative to the original position
    -m : (--max_distance)       The maximum number of nucleotides from the original positions where the new random positions will be generated [Default: 10]
    -n : (--min_distance)       The minimum number of nucleotides from the original positions where the new random positions will be generated [Default: 1]
    -i : (--new_ids)            Column number containing the IDs. When this option is > 0, it will generate new IDs with random numbers. [Default: 0]
    -p : (--id_prefix)          Prefix for the new IDs. the new IDs will contain the prefix followed by a random number. E.g. Sample001 [Default: Sample]


Example:
```unix
## BED file
Rscript R-scripts/Generate_random_coordinates.R \
  -o BED_example                                \
  -t Coordinates.bed                            \
  -f bed                                        \
  -s both                                       \
  -m 15                                         \
  -n 2

## VCF file
Rscript R-scripts/Generate_random_coordinates.R \
  -o VCF_example                                \
  -t Coordinates.vcf                            \
  -f vcf                                        \
  -s both                                       \
  -m 10                                         \
  -n 2                                          \
  -i 5                                          \
  -p DONOR
```

**BED example**: the new coordinates will be generated within a window of 2 to 15 nucleotides from the original coordinates.

**VCF example**: the new coordinates will be generated within a window of 2 to 10 nucleotides from the original coordinates.
Note that the new coordinates will not overlap the original ones. 
In addition, the 5th column contains the IDs, and will be replaced using the prefix 'DONOR', e.g., DONOR001.


_____


### Retrieve motif information from JASPAR 2020

Requires:

  - data.table (CRAN)
  - dplyr (CRAN)
  - jsonlite (CRAN)  
  - purrr (CRAN)  


This scripts connects to an API to retrieve information (TF name, Class, Family, 
References) of all motifs in [JASPAR 2020](http://jaspar.genereg.net/) from a given
taxon and returns a tab-delimited table with this information.

Taxa available:

  - Vertebrates
  - Plants
  - Fungi
  - Insects
  - Nematodes
  - Urochordata

Parameters:

    -o : (--output_directory)  Output folder
    -t : (--taxon)             Taxon available in Jaspar


Example:
```unix
Rscript R-scripts/Retrieve_matrix_information_from_JASPAR2020.R   \
  -o Jaspar_table                                                 \
  -t Vertebrates
```

This script returns a tab-delimited table named *Jaspar_2020_info_{taxon}_table.tab*
wit the following fields:

  - name (TF name)
  - matrix_id (motif ID in JASPAR)
  - class (TF class)
  - family (TF family)
  - tax_group (Taxon)
  - species
  - tax_id
  - type (Data type where this motif was found. E.g., ChIP-seq, SELEX, PBM)
  - uniprot_ids
  - pubmed_ids (Reference to supporting evidence)


_____



### Gene Set Enrichment Analysis with enrichR

Given a gene list, calls *enrichR* and launch a functional enrichment analysis 
with a set of databases. In this script we use the following databases: *KEGG_2019_Human*,
*WikiPathways_2019_Human*, *Panther_2016*, and *GO_Biological_Process_2018*. In 
case you want to use other databases, this should be changed in the code.

The complete list of databases is on the [*enrichR* website](https://maayanlab.cloud/Enrichr/#stats).


Requires:

  - data.table
  - dplyr
  - enrichR
  - ggplot2
  - ggrepel
  - htmlwidgets
  - plotly   


Parameters:

    -g : (--geneset_name)           Name of the gene list. [Default: input_list]
    -l : (--gene_list)              A text file with Entrez gene symbols, one per line, NO header (Mandatory)
    -n : (--N_most_enriched_terms)  Top N most enriched terms to display in the barplots. [Default: 20]
    -o : (--output_directory)       Output folder (Mandatory)
    -p : (--Pvalue_enrichment)      P-value to select relevant terms. [Default: 0.001"] 
    -t : (--title)                  Suffix for output files. [Default: enrichR_analysis]


Example:
```unix
Rscript R-scripts/GSEA_enrichR.R           \
  -l examples/data/dysregulated_genes.txt  \
  -o examples/results/enrichR_results      \
  -g diff_expr_genes                       \
  -n 20                                    \
  -p 0.001                                 \
  -t enrichR_example
```

This script returns:

  - A table for each database ranked by the significance (-log10(Pvalue))
  - Barplots with the N most significant terms above the threshold P-value
  - Plots (static and interactive) with terms ranked by significance. [Example](https://jaimicore.github.io/Doc/Rankplot_terms_Example_enrichR_GO_Biological_Process_2018.html)

___


### Genomic fragments without ambiguous 'N' nucleotides.

Given a genome (e.g., hg19, hg38) returns a BED-like file with the fragments without
ambiguous nucleotides 'N', it means, only contiguous fragments with A, C, G, T.

In a nutshell, for each chromosome this script get the positions of the *N* and 
*{A,C,G,T}*, then calls the R base function *findInterval* to return the intervals
of contiguous *{A,C,G,T}*s. Each chromosome is analyzed independently, taking 
advantage of the parallelization libraries.

We strongly recommend to run this script in a cluster.

For the moment, it only accepts *hg19* and *hg38*, and the chromosomes 1-22, X, Y, and M,
but the code is very easy to expand to include other genomes supported by [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html).


Requires:

  - BSgenome   (Bioconductor)
  - doParallel
  - dplyr
  - foreach


Parameters:

    -g : (--genome_version)     Genome version. [Default: hg38]
    -c : (--cores)              Number of cores to parallelize the process. [Default: 3]
    -o : (--output_directory)   Output directory to export the results (Mandatory)


Example:
```unix
Rscript R-scripts/NonAmbiguous_nucleotides.R   \
  -o .                                         \
  -g hg38                                      \
  -c 25
```

This script returns a BED-like table with the following columns:

  - Chromosome
  - Start
  - End
  - Segment width
  - Fraction of the segment with respect to the total number of non-ambiguous nucleotides in the entire chromosome

```unix
Chr   Start   End     Width   Percent
chr1  10001   207666  197665  0.0008576195
chr1  257667  297968  40301   0.0001748561
chr1  347969  535988  188019  0.0008157679
```
___


### Motif friseur

A small script to trim Transcription Factor Position Frequency Matrix

Supported motif formats:

  - jaspar
  
Other motifs formats will be added upon request or if they are required.


Requires:

  - universalmotif  (Bioconductor)


Parameters:

    -m : (--input_motif)       Path to input motif file. (Mandatory)
    -n : (--format)            Input motif format. (Mandatory)
    -o : (--output_directory)  Output directory to export the trimmed motifs. If not indicated, the trimmed motifs are exported in the same folder with the extension *.trimmed*
    
    -b : (--both)              Trim *b* nucleotides in both sides 
    
    -l : (--left)              Trim *l* nucleotides in left side
    -r : (--right)             Trim *r* nucleotides in right side
    
    -f : (--from)              Keep nucleotides starting from *f* position 
    -t : (--to)                Keep nucleotides until *t* position 
    

Example:
```unix
Rscript R-scripts/Motif_Friseur.R   \
  -m examples/data/ZNF506.jaspar    \
  -n jaspar                         \
  -o examples/results/Motif_friseur \
  -f 4                              \
  -t 12
```

Trimming the ZNF506 motif from [JASPAR](http://jaspar.genereg.net/matrix/UN0198.1/),
the trimmed motif contains the columns 4 to 12.



___


### Generate randomized (weighted) networks

Given an input weighted network, returns a randomized version without duplicated
targets. The input network must have at least two columns: 1) Regulator and 2) Targets,
optionally it may  contain a third column with Weights.

Randomization approaches:

  - *Simple*: the targets in the network are pooled and for each regulator we sample X targets, where X is the original number of targets.
  - *List*: similar to the simple approach, but using an user-provided list of targets.
  - *Shuffle*: re-distribute the original targets avoiding duplicated entries on each regulator. This approach preserves the number of times each target gene appears in the input network.


In the *shuffle* mode, if a gene (e.g., *TP53*) is regulated by 15 genes in the input
network, the resulting randomized network will contain *TP53* regulated by 15 
distinct genes, this is important to preserve the original connectivity of the
input network but we may lost the biological relevance of the connections.


Requires:

  - data.table   
  - dplyr
  - purrr


Example of input network file:

```unix
Gene    Target    Weight
EP300	GNA15	  0.678
EP300	PIK3CG	  0.899
EP300	EPO	  0.882
SRY	AP1G2	  0.663
GATA3	NR5A1	  0.505
GATA3	CENPH	  0.63
```


Example of input target list file:

```unix
GNA15
PIK3CG
EPO
AP1G2
NR5A1
CENPH
```

Parameters:

    -n : (--network_file)	Network file. (Mandatory)
    -o : (--output_directory)		Output directory to export the results (Mandatory)
    -s : (--suffix)		Suffix added in the output file name. [Default: Random_network_1]
    -m : (--mode)		Indicates how the random network should be created. [Default: simple] [Options: simple, list, shuffle] 
    -w : (--permute_weights)		Indicates if the weights should be permuted or no.. [Default: 0] [Options: 0 | 1] 
    -l : (--list)		A text file containing the names (a name per line, no header) of all the targets that will be sampled to create the random network.


Example:
```unix

## Simple mode
Rscript R-scripts/Randomize_weighted_network.R  \
  -n examples/data/Weighted_net_example.txt     \
  -o examples/results/Random_networks           \
  -s Xseq_random_net_1                          \
  -m simple
  
  
## List mode
Rscript R-scripts/Randomize_weighted_network.R  \
  -n examples/data/Weighted_net_example.txt     \
  -o examples/results/Random_networks           \
  -s Xseq_random_net_1                          \
  -m list                                       \
  -l examples/data/Target_list_example.txt 
  
  
## Shuffle mode
Rscript R-scripts/Randomize_weighted_network.R  \
  -n examples/data/Weighted_net_example.txt     \
  -o examples/results/Random_networks           \
  -s Xseq_random_net_1                          \
  -m shuffle
```


This script returns a tab file containing a randomized version of the input
network, same structure as the input network file (see above).

___


### Shuffle table columns and/or rows 

Given an input table file, shuffle the indicated column(s) and/or row(s).

Requires:

  - data.table


Parameters:

    -c : (--columns)            Indicates the columns to shuffle. Example: -c 1,2. If this value is 0, shuffle all columns, if this is not indicated, none column is shuffled.
    -r : (--rows)               Indicates the rows to shuffle. Example: -c 1,2. If this value is 0, shuffle all rows, if this is not indicated, none row is shuffled.
    -h : (--header)             Indicates wheter the input table has a header [Y | N]. [Default: Y]
    -o : (--output_directory)   Output directory to export the results (Mandatory)
    -i : (--input_table)        Input table
    -p : (--prefix)             A prefix added to the output file name. [Default: Shuffled]


Example:
```unix
Rscript R-scripts/Shuffle_table.R   \
  -i examples/data/TFBSs.bed        \
  -o examples/results/Shuffle_table \
  -p Shuffled_coordinates           \
  -c 1,4,5
```


This script returns a table with the shuffled columns/rows.

___


### Group UniBind datasets given a keyword

Concatenates all [*UniBind*](https://unibind.uio.no/) datasets associated to a given keyword (e.g., lung, cancer, ESC). Returns a BED file with all the associated TFBSs. 

First, download the complete UniBind results (given the genome and collection),
if this file exist in the output directory, then this step is skipped. Note that
downloading the results make take several minutes, this depends on your internet
collection.

The keyword is searched on the column *title* in the UniBind metadata table.

Please contact the UniBind developers to have access to the metadata information: r.r.puig AT ncmm.uio.no

Currently this scripts support only one genome (human), however, UniBind has data
of 9 species, the remaining ones will be added soon.

Requires:

  - data.table
  - dplyr
  - purrr
  

Parameters:

    -w : (--word)               Keyword used to return dataset containing such word in their description (Mandatory)
    -g : (--genome)             Available genomes in UniBind: human.
    -c : (--collection)         Dataset quality: robust | permissive. [Default: Robust]
    -o : (--output_directory)   Output directory to export the results (Mandatory)
    -a : (--annotation_table)   UniBind metadata table (Mandatory)
    -p : (--prefix)             A prefix added to the output file name. [Default: Shuffled]


Example:
```unix
Rscript R-scripts/Shuffle_table.R   \
  -i examples/data/TFBSs.bed        \
  -o examples/results/Shuffle_table \
  -p Shuffled_coordinates           \
  -c 1,4,5
```

Output:

This script returns a BED file containing the datasets that matched the keyword in
their description.
___
