#####################
## Load R packages ##
#####################
message("; Loading required packages")
required.libraries <- c("data.table",
                        "dplyr",
                        "ggplot2",
                        "optparse",
                        "purrr",
                        "tidyr"
                        )


for (lib in required.libraries) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}



####################
## Read arguments ##
####################

## 02-08-2021 This function was taken and adapted from from :
## https://github.com/pablo-gar/Rmods/blob/master/optparse_template.R

getOpts <- function(positional_arguments=F) {
  
  ## Specify the  type of options, this will aid for checking the type of opts
  required  <- c('input', 'output') ## Will  throw an error if these are not  given
  files     <- c('input')           ## Will throw an error if the path to this options does not exist
  out_files <- c('output')          ## Will throw an error if the parent dir to these files don't exist
  dirs      <- c()                  ## Will check that these folders exist
  
  
  ## Create arg list
  option_list <- list(
    
    make_option(c('-o', '--outoput'), type = 'character', help = 'Path to output file'),
    
    make_option(c('-i', '--input'), type = 'character', help = 'Path to input file'),
    
    make_option(c("-g", "--dysregulated_gene_threshold"), type="numeric", default = 0.5, 
                help = "Threshold to select dysregulated genes. [Default \"%default\"] ", metavar = "number"),
    
    make_option(c("-d", "--output_directory"), type = "character", default = NULL, 
                help = "Output directory to export the results (Mandatory)", metavar = "character")
  )
  
  
  ############################################################
  ## DO NOT MODIFY AFTER THIS UNTIL THE END OF THE FUNCTION ##
  ############################################################
  opt_parser <- OptionParser(usage       = 'usage: %prog [options]', 
                             option_list = option_list, 
                             description = 'Calculate median gene expression across regions')
  
  opt <- parse_args(opt_parser, positional_arguments=positional_arguments)
  
  if (positional_arguments) {
    opt_check <- opt$options
  } else {
    opt_check <- opt
  }
  
  
  ## Checking for essential arguments ##
  for (i in required) {
    if (is.null(opt_check[i][[1]])) {
      stop('"--', i, '" is a required argument, run with "-h" for help')
    }
  }
  
  
  ## Checking files exist
  for (i in files) {
    if (!file.exists(opt_check[i][[1]])) {
      stop('"--', i, '" "', opt_check[i][[1]], '" file does not exist')
    }
  }
  
  
  ## Checking that we can write out files
  for (i in out_files) {
    if (!dir.exists(dirname(opt_check[i][[1]]))) {
      stop('"--', i, '" "', opt_check[i][[1]], '" parent folder does not exist')
    }
  }
  
  ## Checking  dirs exists()
  for (i in dirs) {
    if (!dir.exists(opt_check[i][[1]])) {
      stop('"--', i, '" "', opt_check[i][[1]], '" folder does not exist')
    }
  }
  
  return(opt)
  
}

## Retrieve and Check the arguments
message("; Reading arguments from command-line")
opts <- getOpts(positional_arguments = F)



##############################
## Set input variable names ##
##############################
results.dir <- opt$output_directory

## Numeric
dysreg.gene.thr <- as.numeric(opt$dysregulated_gene_threshold)

## Default
min.nb.signif.samples <- 2
do.data.type <- list()



###########################
## Create results folder ##
###########################
message("; Creating output folders")
out.folders <- list()
out.folders[["plots"]]  <- file.path(results.dir, "plots")
out.folders[["tables"]] <- file.path(results.dir, "tables")
out.folders[["RData"]]  <- file.path(results.dir, "Rdata")
no.message <- sapply(out.folders, dir.create, recursive = TRUE, showWarnings = FALSE)



##########################
## Add sections here ...