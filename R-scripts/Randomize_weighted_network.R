#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "doParallel",
                        "dplyr",
                        "foreach",
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
              help = "Output directory to export the tables (Mandatory)", metavar = "character"),
  
  make_option(c("-n", "--network_file"), type = "character", default = "NULL", 
              help = "(Mandatory input file) File containing the network, must include at least 2 columns: 1) Gene, 2) Target, and optionally 3) Weight ", metavar = "character"),  
  
  make_option(c("-r", "--random_networks"), type = "numeric", default = 10, 
              help = "Number of random networks. [Default \"%default\"] ", metavar = "number"), 
  
  make_option(c("-c", "--cores"), type = "numeric", default = 2, 
              help = "Number of cores to parallelize. [Default \"%default\"] ", metavar = "number"), 
  
  );
message("; Reading arguments from command line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


########################
## Set variable names ##
########################
results.dir <- opt$output_directory
net.file    <- opt$network_file
nb.rand.net <- as.numeric(opt$random_networks)
nb.cores    <- as.numeric(opt$cores)


###########
## Debug ##
###########
net.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/examples/data/Weighted_net_example.txt"
nb.rand.net <- 2
nb.cores <- 3


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument (Output directory): results.dir ")
  
} else if (!exists("net.file")) {
  stop("Missing mandatory argument: net.file ")
  
}
dir.create(file.path(results.dir, "Random_networks"), showWarnings = F, recursive = T)


#######################
## Read network file ##
#######################
net <- fread(net.file)

## Check if the network contains weights or not
w.flag <- 0
if (ncol(net) >= 3) {
  w.flag <- 1
  colnames(net)[3] <- "Weight"
}

colnames(net)[1] <- "Gene"
colnames(net)[2] <- "Target"


## Network represented as a list, where each element corresponds to
## a dataframe of the 'Gene' columns
net.list  <- split(net, f = net$Gene)

## Sort the networks by decreasing size
# gen.net.size.order <- order(unlist(lapply(net.list, nrow)), decreasing = T)
# net.list           <- net.list[gen.net.size.order]

registerDoParallel(nb.cores)
new.targets.list <- foreach(i = 1:nb.rand.net) %dopar% {
  
  ## The universe of target genes (many of them are repeated because are target of many genes)
  ## We keep this number to maintain exactly the same number of each target in the random network
  all.targets <- net$Target
  
  ## Iterate over each Gene in the network
  no.dup.tx.flag <- 0
  while (!no.dup.tx.flag) {
    
    new.targets <- lapply(net.list, function(l){
      
      message("; Number of targets to reallocate: ", length(all.targets))
      
      ## Network size/Number of targets
      nb.entries <- nrow(l)
      
      ## Shuffle the names on each iteration
      all.targets <<- sample(all.targets)
      
      ## Get a new set of non-duplicated target genes
      new.genes <- unique(all.targets)[1:nb.entries]
      
      
      new.genes.ind <- match(new.genes, all.targets)
      all.targets <<- all.targets[-new.genes.ind]
      
      new.genes
    })
    
    dup.sum <- sum(unlist(lapply(lapply(new.targets, duplicated), sum)))
    
    message("; Duplicated entries in the new random network: ", dup.sum)
    no.dup.tx.flag <- ifelse(dup.sum == 0, yes = 1, no = 0)
  }
  new.targets
}



new.targets.df <- lapply(lapply(new.targets.list, lapply, data.frame), purrr::map_dfr, data.table)

lapply(new.targets.df, function(l){
  data.frame(Gene   = net$Gene,
             Target = l,
             Weight = sample(net$Weight))
})



# colnames(new.targets.df) <- "Target"
# 
# new.targets.df <- purrr::map(new.targets.list, data.frame)

