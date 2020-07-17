#####################
## Load R packages ##
#####################
required.libraries <- c("biomaRt",
                        "optparse")
for (lib in required.libraries) {
  if (!require(lib, character.only=TRUE)) {
    install.packages(lib)
    library(lib, character.only=TRUE)
  }
}
## Status: revised (17-02-2020) 


####################
## Read arguments ##
####################
option_list = list(
  make_option(c("-g", "--human_genome_version"), type="character", default=NULL, 
              help="Human genome version. Options: hg19 | hg38", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default=NULL, 
              help="Output directory to export the tables", metavar="character")
  
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

results.dir <- opt$output_directory
hg.version <- opt$human_genome_version

## Numeric version
hg.version <- as.numeric(gsub(hg.version, pattern = "\\D+", replacement = ""))

## Rename genome version
if(hg.version == "hg19"){
  hg.version <- 37
} else if(hg.version == "hg38"){
  hg.version <- 38
}


## ================================= ##
## Get gene annotations from ENSEMBL ##
## ================================= ##
message("; Retrieving gene annotations from ENSEMBL")

## NOTE: the RNA-seq expression contains gene and transcript IDs
## The transcript IDs should be mapped to their corresponding genes using the ENSEMBL biomart functions
## The table must be reduced to included protein-coding genes only
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = hg.version) # GRCh=37

## Get the gene and transcript IDs for all human genes
## Filter the selection by protein coding entries
## Genome version: hg38
all.genes <- getBM(attributes=c('ensembl_gene_id',
                                'ensembl_transcript_id',
                                'ensembl_peptide_id',
                                'hgnc_symbol',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                "gene_biotype",
                                'entrezgene_id'),
                   filters = c("biotype"),
                   values = list("protein_coding"),
                   mart = ensembl)
colnames(all.genes)[1:3] <- c('gene_id', 'transcript_id', 'peptide_id')
# View(all.genes)


## =============================== ##
## Filter the table by chromosomes ##
## =============================== ##

## Keep genes in chr 1-22, X, Y, MT
## Other chromosomes (e.g., CHR_HSCHR6_MHC_MANN_CTG1) are removed.
## The names refer to genomic sequence that differs from the genomic DNA on the primary assembly.
## These alternate sequences come in two types: Allelic sequence (haplotypes and novel patches) and fix patches.
## See: http://www.ensembl.info/2011/05/20/accessing-non-reference-sequences-in-human/
selected.chr <- c(1:22, "MT", "X", "Y")
selected.chr <- as.character(selected.chr)
all.genes <- all.genes[all.genes$chromosome_name %in% selected.chr, ]


## ======================== ##
## Update the ENSEMBL table ##
## ======================== ##

## Add the gene name of the following ENSEMBL IDs. Information retrieved from genecards.
## According to targetScan, these genes are miRNA targets, however they didn't match with the ENSEMBL table. In the table, they have not gene names, but they are considered as protein coding genes.
## The names, obtained from genecards, will be inserted in the table.

genes.to.add <- c("SGK494")
ids.to.add <- c("ENSG00000167524")
names(ids.to.add) <- genes.to.add

# http://www.genecards.org/cgi-bin/carddisp.pl?gene=PLPPR1&keywords=LPPR1
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=PLPPR4&keywords=LPPR4
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=PLPPR5&keywords=LPPR5
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=LZTS3&keywords=LZTS3
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=MAP3K21&keywords=MLK4
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=MAP3K20&keywords=MLTK
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=STK26&keywords=MST4
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=NSG1&keywords=nsg1
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=NSG2&keywords=nsg2
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=ACP7&keywords=PAPL
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=SELENOK&keywords=SELK
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=SELENOM&keywords=SELm
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=SELENOT&keywords=selt
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=SELENOV&keywords=SELV
# http://www.genecards.org/cgi-bin/carddisp.pl?gene=SGK494&keywords=SGK494

for(i in 1:length(ids.to.add)){
  all.genes[which(all.genes$gene_id == ids.to.add[i]) ,]$hgnc_symbol <- names(ids.to.add[i])
}


## Remove entries without gene name and entrezgene ID
## These genes cannot be mapped, therefore will not be considered in the further analysis.
# all.genes2 <- all.genes
# all.genes <- all.genes2
dim(all.genes)
all.genes <- all.genes[all.genes$hgnc_symbol != "",]
all.genes <- all.genes[all.genes$peptide_id != "",]
all.genes <- all.genes[!is.na(all.genes$entrezgene),]
all.genes <- unique.data.frame(all.genes)
dim(all.genes)

# gene.id.tab <- unique(all.genes[,c(1,4)])


## ================================================== ##
## Export ENSEMBL annotations in tab and Rdata format ##
## ================================================== ##
message("; Exporting the protein coding genes from the ENSEMBL table")

## Remove the transcript_id column
ENSEMBL.protein.coding.genes <- all.genes#[,c(-2)]
# ENSEMBL.protein.coding.genes <- unique.data.frame(ENSEMBL.protein.coding.genes)

ensemble.table.file <- file.path(results.dir, "ENSEMBL_protein_coding_genes.tab")
ensemble.table.file.rdata <- file.path(results.dir, "ENSEMBL_protein_coding_genes.Rdata")

write.table(ENSEMBL.protein.coding.genes, file = ensemble.table.file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
save(ENSEMBL.protein.coding.genes, file = ensemble.table.file.rdata)

