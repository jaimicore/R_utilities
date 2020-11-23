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



####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-o", "--output_directory"), type = "character", default = NULL, 
              help = "Output directory to export the table (Mandatory)", metavar = "character"),
  
  make_option(c("-i", "--input_table"), type = "character", default = NULL, 
              help = "Input table. (Mandatory) ", metavar = "character"),  
  
  make_option(c("-h", "--header"), type = "character", default = "Y", 
              help = "Indicates wheter the input table has a header [Y | N]. [Default \"%default\"]. ", metavar = "character"),  
  
  make_option(c("-c", "--columns"), type = "character", default = NULL, 
              help = "Indicates the columns to shuffle. Example: -c 1,2. If this value is 0, shuffle all columns, if this is not indicated, none column is shuffled.", metavar = "character"),  
  
  make_option(c("-r", "--rows"), type = "character", default = NULL, 
              help = "Indicates the rows to shuffle. Example: -c 1,2. If this value is 0, shuffle all rows, if this is not indicated, none row is shuffled.", metavar = "character"),
  
  make_option(c("-p", "--prefix"), type = "character", default = "Shuffled", 
              help = "A prefix added to the output file name. [Default \"%default\"].", metavar = "character")
);
message("; Reading arguments from command line")
opt_parser = OptionParser(option_list = option_list, add_help_option = F);
opt = parse_args(opt_parser);


########################
## Set variable names ##
########################
message("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
results.dir    <- opt$output_directory
input.tab.file <- opt$input_table
header.flag    <- tolower(opt$header)
columns        <- opt$columns
rows           <- opt$rows
prefix         <- opt$prefix


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument (Output directory): results.dir ")
  
} else if (!exists("input.tab.file")) {
  stop("Missing mandatory argument: input.tab.file ")
  
}
dir.create(results.dir, showWarnings = F, recursive = T)


###########
## Debug ##
###########
# input.tab.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/examples/data/TFBSs.bed"
# header.flag    <- "y"
# prefix         <- "Example"
# rows           <- NULL
# columns        <- "1,2"
# 
# columns        <- "1"
# rows           <- "2"
# 
# columns        <- "0"
# rows           <- "0"


######################
## Read input table ##
######################

## Parse header parameters
if (header.flag == "y") {
  header.flag <- T
} else if (header.flag == "n") {
  header.flag <- F
} else {
  stop("; Header must be indicated as Y|N")
}

message("; Reading input table: ", input.tab.file)
input.tab <- fread(input.tab.file, header = header.flag)
input.tab <- data.frame(input.tab)

#################################
## Shuffle columns or/and rows ##
#################################

## Columns
columns.flag <- 0
if (!is.null(columns)) {
  
  columns <- as.numeric(unlist(strsplit(columns, split = ",")))
  columns.flag <- 1
  
  ## When the value is 0, all columns will be shuffled
  if (columns == 0) {
    columns <- 1:ncol(input.tab)
  }
}

## Rows
rows.flag <- 0
if (!is.null(rows)) {
  
  rows <- as.numeric(unlist(strsplit(rows, split = ",")))
  rows.flag <- 1
  
  ## When the value is 0, all rows will be shuffled
  if (rows == 0) {
    rows <- 1:nrow(input.tab)
  }
}


## Shuffle columns when indicated
if (columns.flag) {
  
  message("; Shuffling columns")
  for (i in columns) {
    input.tab[, i] <- sample(as.vector(unlist(input.tab[, i])))
  }
}


## Shuffle rows when indicated
if (rows.flag) {
  
  message("; Shuffling rows")
  for (i in rows) {
    
    new.row           <- data.frame(sample(input.tab[i,]))
    colnames(new.row) <- colnames(input.tab)
    input.tab[i,]     <- new.row
  }
}


###########################
## Export shuffled table ##
###########################
output.tab <- data.table(input.tab)
output.tab.file <- file.path(results.dir, paste0(prefix, "_", basename(input.tab.file)))
message("; Reading input table: ", output.tab.file)
fwrite(output.tab, sep = "\t", row.names = F, col.names = header.flag)