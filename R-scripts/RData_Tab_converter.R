#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "optparse")
for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


## How to run: 
## 
## From Table to Rdata: Rscript RData_Tab_converter.R -i TFBSs_Unibind_EZ_hg19.RData -from RData -to tab -o . -e bed
##
## From Rdata to Table: Rscript RData_Tab_converter.R -i TFBSs_Unibind_EZ_hg19.bed -from tab -to RData -o .


####################
## Read arguments ##
####################
option_list = list(
  make_option( c("-f", "--from"), type = "character", default = NULL, help = "Input data type: tab | RData", metavar = "character"),
  make_option( c("-t", "--to"), type = "character", default = NULL, help = "Output data type: tab | RData", metavar = "character"),
  make_option( c("-i", "--input_file"), type = "character", default = NULL, help = "Input file", metavar = "character"),
  make_option( c("-o", "--out_folder"), type = "character", default = ".", help = "Output folder", metavar = "character"),
  make_option( c("-r", "--rownames"), type = "numeric", default = 0, help = "Print rownames in table: 0|1. Default: 0", metavar = "number"),
  make_option( c("-l", "--colnames"), type = "numeric", default = 0, help = "Print colnames in table: 0|1. Default: 0", metavar = "number"),
  make_option( c("-e", "--tab_extension"), type = "character", default = "tab", help = "Output table extension. Example: BED, tab, txt", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

input.file    <- opt$input_file
input.format  <- opt$from
output.format <- opt$to
out.folder    <- opt$out_folder
rownames.df   <- opt$rownames
colnames.df   <- opt$colnames
tab.ext       <- opt$tab_extension

dir.create(out.folder, showWarnings = F, recursive = T)
message("; Input table: ", input.file)

## From table to RData
output.file <- NULL
if (input.format == "tab" & output.format == "RData") {
  
  ## Read input table
  input.table <- fread(input.file)
  
  ## Export the object as RData 
  output.file <- gsub(x = basename(input.file), pattern = "\\..+$", replacement = "\\.RData")
  output.file <- file.path(out.folder, output.file)
  save(input.table, file = output.file)
  
  ## From RData to table
} else if (input.format == "RData" & output.format == "tab") {
  
  ## Load Rdata object
  rdata.file <- load(input.file)
  df <- get(rdata.file)
  
  ## Export table
  output.file <- gsub(x = basename(input.file), pattern = "\\.RData$", replacement = paste0("\\.",tab.ext))
  output.file <- file.path(out.folder, output.file)
  fwrite(df, file = output.file, sep = "\t", col.names = as.logical(colnames.df), row.names = as.logical(rownames.df))
  
} else {
  
  stop("; Input and output format not recongized: [tab | RData]")
  
}

message("; Conversion ready: ", output.file)
