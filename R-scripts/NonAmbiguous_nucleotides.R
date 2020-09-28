#####################
## Load R packages ##
#####################
required.libraries <- c("BSgenome",
                        "doParallel",
                        "dplyr",
                        "foreach",
                        "optparse")

for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-o", "--output_directory"), type = "character", default = NULL, 
              help = "Output directory to export the results (Mandatory)", metavar = "character"),
  
  make_option(c("-g", "--genome_version"), type = "character", default = "hg38", 
              help = "Human genome version. [hg38|hg19 : Default \"%default\"] ", metavar = "character"),
  
  make_option(c("-c", "--cores"), type = "numeric", default = 3, 
              help = "Number of cores to paralellize. [Default \"%default\"] ", metavar = "number")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


######################
## Define functions ##
######################
find.non.ambiguous.regions <- function(Seqname = "Chr1",
                                       DNAseq = "NNNTACGATCTACGCGATTATGCCGATTAAGCACGTACNCTACGATCGATCGATCGATGATCGCNNGCATCTGATCGTACATGCTGAGCACGTAGAGCNN") {
  
  DNAseq <- toupper(DNAseq)
  
  ## Separate the string into a vector of N characters
  DNAseq.split <- unlist(strsplit(DNAseq, split = ""))
  
  ## Get positions with {A,C,G,T} and {N}
  acgt.nt.pos <- which(DNAseq.split != "N")
  n.nt.pos    <- which(DNAseq.split == "N")
  
  acgt.length <- length(acgt.nt.pos)
  
  ## Find the intervals
  ##
  ## Get the intervals between a pair of increasing numeric sequences
  ## in this case, we prove the vector with the positions of {A,C,G,T} and compare it
  ## againts the N's positions.
  ## Then, simplify the results to return a BED-like format
  interval.gaps <- data.frame(cbind(acgt.nt.pos,
                                    findInterval(acgt.nt.pos, n.nt.pos))) %>% 
    dplyr::rename(Position = "acgt.nt.pos",
                  Scaffold = "V2") %>% 
    group_by(Scaffold) %>% 
    summarise(Chr   = Seqname,
              Start = min(Position),
              End   = max(Position)) %>%
    mutate(Width = End - Start) %>%
    mutate(Percent = round(Width/acgt.length, digits = 10)) %>% 
    select(Chr, Start, End, Width, Percent)
  
  return(interval.gaps)
}


########################
## Set variable names ##
########################
results.dir    <- opt$output_directory
genome.version <- opt$genome_version
ncores         <- opt$cores


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument: results.dir ")
  
}

## Debug
# genome.version <- "hg38"


###########################
## Select genome version ##
###########################
if (genome.version == "hg38") {
  gv <- "BSgenome.Hsapiens.UCSC.hg38"
  
} else if (genome.version == "hg19") {
  gv <- "BSgenome.Hsapiens.UCSC.hg19"

} else {
  stop("; Human genome version non recognized. Options: hg19, hg38")
}
## Load the selected genome
library(gv, character.only = T)


## Install the selected genome
# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install(gv)
# # "BSgenome.Hsapiens.UCSC.hg19" 
# # "BSgenome.Hsapiens.UCSC.hg38"


#######################################################
## Get the regions without ambiguous 'N' nucleotides ##
#######################################################
chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
# chromosomes <- paste0("chr", c(22, "Y"))
# chr <- "chr22"


registerDoParallel(ncores)
non.abiguous.regions.df <- NULL
non.abiguous.regions.df <- foreach(i = seq_along(chromosomes), .combine = rbind) %dopar% {

  chr <- chromosomes[i]
  message("; Processing chromosome ", chr)
  
  chrom.seq <- as.character(BSgenome.Hsapiens.UCSC.hg38[[chr]])
  
  find.non.ambiguous.regions(Seqname = chr, DNAseq  = chrom.seq)
  
}


## Export table with non-ambiguous regions
non.abiguous.regions.df.file <- file.path(results.dir, paste0("Non_ambiguous_regions_", genome.version,".bed"))
write.table(non.abiguous.regions.df, file = non.abiguous.regions.df.file, quote = F, row.names = F, col.names = T)
