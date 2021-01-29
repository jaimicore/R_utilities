#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "optparse",
                        "purrr")

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
              help = "Output directory to export the table (Mandatory)", metavar = "character"),

  make_option(c("-w", "--word"), type = "character", default = "breast", 
              help = "Keyword used to return dataset containing such word in their description (Mandatory)", metavar = "character"),
  
  make_option(c("-g", "--genome"), type = "character", default = "human", 
              help = "Available genomes: human", metavar = "character"),
  
  make_option(c("-c", "--collection"), type = "character", default = "Robust", 
              help = "Dataset quality: robust | permissive", metavar = "character"),
  
  make_option(c("-a", "--annotation_table"), type = "character", default = NULL, 
              help = "UniBind metadata table (Mandatory)", metavar = "character")
);
message("; Reading arguments from command line")
opt_parser = OptionParser(option_list = option_list, add_help_option = F);
opt = parse_args(opt_parser);


## Upper case for first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


########################
## Set variable names ##
########################
results.dir       <- opt$output_directory
keyword           <- opt$word
genome            <- firstup(opt$genome)
collection        <- firstup(opt$collection)
metadata.tab.file <- opt$annotation_table



###########
## Debug ##
###########
# results.dir     <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/Test1"
# keyword         <- "breast"
# genome          <- "Human"
# collection      <- "Robust"
# metadata.tab.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/Test1/UniBind_2021/Homo_sapiens/hg38/metadata.tsv"


#####################
## Assign assembly ##
##                 ##  
## To be completed ##
##                 ##
#####################
if (genome == "Human") {
  
  genome   <- "Homo_sapiens"
  assembly <- "hg38"
  
} else {
  stop("Genome not found: ", genome, ". Supported ones: human")
}


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument (Output directory): results.dir ")
  
} else if (!exists("keyword")) {
  stop("Missing mandatory argument: keyword ")
  
}


######################
## Global variables ##
######################
unibind.datasets.dir   <- file.path(results.dir, "UniBind_2021", genome, assembly)
unibind.download.link  <- paste0("https://unibind.uio.no/static/data/20201216/bulk_", collection, "/", genome, "/damo_", assembly, "_TFBS_per_TF.tar.gz")
unibind.datasets.gz    <- file.path(results.dir, "UniBind_2021", genome, assembly, paste0("damo_", assembly, "_TFBS_per_TF.tar.gz"))
unibind.selection.dir  <- file.path(results.dir, "UniBind_2021", "Subset_by_keyword", keyword)

## Create output directories
message("; Creating output directories")
dir.create(unibind.datasets.dir, showWarnings = F, recursive = T)
dir.create(unibind.selection.dir, showWarnings = F, recursive = T)


###############################
## Download UniBind datasets ##
###############################
if (!file.exists(unibind.datasets.gz)) {
  message("; Downloading UniBind data")
  unibind.temp <- tempfile(fileext = '.tar.gz')
  download.file(url = unibind.download.link, destfile = unibind.temp, quiet = FALSE, method = 'curl')
  message("; Download complete")
  
  ## Move to copy to the UniBind directory
  file.copy(from = unibind.temp, to = unibind.datasets.gz)
  file.remove(unibind.temp)
  
  ## Decompress the file
  untar(unibind.datasets.gz, compressed = 'gzip', exdir = unibind.datasets.dir)
}


#########################
## Read metadata table ##
#########################
message("; Reading metadata table")
species.name <- gsub(genome, pattern = "_", replacement = " ")
metadata.tab <- fread(metadata.tab.file) %>% 
                  dplyr::filter(species == species.name)

## Select entries matching the keyword
message("; Finding datasets matching the keyword: ", keyword)
metadata.selection <- metadata.tab[grepl(metadata.tab$title, pattern = keyword, ignore.case = T), ]
metadata.selection <- metadata.selection[!duplicated(metadata.selection$UBID), ]

## Check that the selection is not an empty table
if (nrow(metadata.selection) == 0) {
  stop("; The keyword ", keyword, " was not found in the metadata table.")
}


######################################
## Generate table with all datasets ##
######################################
message("; Obtaining the list of all the available BED files")
## Get the name of all the BED files in their respective TF folders
TFs.folder   <- file.path(unibind.datasets.dir, "TFBS_per_TF")
TF.files.tab <- sapply(list.files(TFs.folder), function(TF){
  
  list.files(file.path(TFs.folder, TF))
  
})


## Create a dataframe with all the BED files and their associated TFs
TFBS.bed.files.tab <- data.frame(TF   = rep(names(TF.files.tab), sapply(TF.files.tab, length)),
                                 File = unlist(TF.files.tab)) %>% 
                        mutate(TFBS_BED_file = file.path(file.path(TFs.folder, TF, File)))

# TFBS.bed.files.tab <- data.frame(TF   = rep(names(TF.files.tab), sapply(TF.files.tab, length)),
#                                  File = unlist(TF.files.tab),
#                                  UBID = gsub(unlist(TF.files.tab), pattern = "_MACS_.+$", replacement = "", ignore.case = T, perl = T)) %>% 
#   mutate(TFBS_BED_file = file.path(file.path(TFs.folder, TF, File)))


##############################
## Select relevant datasets ##
##############################
message("; Selecting relevant datasets")

## First match the UniBind ID (from the metadata table) in the BED file names
selected.bed.files <- lapply(metadata.selection$UBID, function(ub){
                        TFBS.bed.files.tab$File[grep(x = TFBS.bed.files.tab$File, pattern = ub)]
                      })
selected.bed.files <- unlist(selected.bed.files)

## Keep the names that matched the UBID 
## NOTE: the metadata contains both permissive and robust datasets, therefore many 
##       datasets intially selected using the keyword, may not be in the robust
##       collection
TFBS.bed.files.tab <- TFBS.bed.files.tab[TFBS.bed.files.tab$File %in% selected.bed.files, ]


## Concat all the files in a single table
message("; Concatenating all relavant datasets in a single table")
metadata.selection
TFBS.bed.df <- TFBS.bed.files.tab$TFBS_BED_file %>% 
                map_df(~fread(.))



########################################################################################
## Export the concatenated set of TFBSs whose dataset description contain the keyword ##
########################################################################################
message("; Exporting concatenated set of TFBSs as a BED file")
TFBS.bed.file <- file.path(unibind.selection.dir, paste0("UniBind_TFBS_selection_keyword_", keyword, "_", genome, "_", assembly, ".bed"))
fwrite(TFBS.bed.df, file = TFBS.bed.file, sep = "\t", row.names = F, col.names = F)