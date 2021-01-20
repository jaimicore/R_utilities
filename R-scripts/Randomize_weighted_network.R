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
  
  make_option(c("-c", "--cores"), type = "numeric", default = 2, 
              help = "Number of cores to parallelize. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-e", "--exclusive"), type = "numeric", default = FALSE, 
              help = "Indicates whether the random targets must not being the original ones. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-s", "--suffix"), type = "character", default = NULL, 
              help = "Suffix name for the random network file. [Default \"%default\"]. By default the output files follow this format: Random_network_{it}.tab, where {it} correspond to the number of random network. When the suffix option is indicated, change {it} will take the suffix value. Note that the suffix is applied to all file, so this option is recommendended when only one random net work is generated.", metavar = "character")
  
  );
message("; Reading arguments from command line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


########################
## Set variable names ##
########################
results.dir     <- opt$output_directory
net.file        <- opt$network_file
exclude.targets <- opt$exclusive
nb.rand.net     <- as.numeric(opt$random_networks)
nb.cores        <- as.numeric(opt$cores)
suffix.net.name <- opt$suffix

nb.rand.net     <- 1
nb.cores        <- 1

## Number of cores must be at most the number of generated random networks
nb.cores <- ifelse(nb.cores > nb.rand.net, yes = nb.rand.net, no = nb.cores)


########################################################
## Find a Gene-Target association to exchange targets ##
########################################################
# net.df = rand.net.df
# gene   = "hsa-mir-1234-3p"
# target = "LPP"
# no.dup = T
# no.ori = T
find.new.gene.target <- function(net.df = NULL,
                                 gene   = NULL,
                                 target = NULL,
                                 no.dup = T,
                                 no.ori = F) {
  
  gene.tx.subset   <- subset(net.df, Gene == gene & Target == target)
  gene.subset      <- subset(net.df, Gene == gene)
  
  nb.g.tx.entries  <- nrow(gene.tx.subset)
  genes.set        <- unique(net.df$Gene)
  half.set.size    <- floor(length(genes.set)/2)
  problematic.flag <- 0
  tx.in.net        <- gene.subset$Target
  ori.tx.in.net    <- gene.subset$Ori_tx
  from.entry       <- gene.tx.subset$ID[1]
  
  ## Check if the Gene-Target association exists in the gene network
  if (nb.g.tx.entries) {
    
    ## Constraint combinations
    if (no.dup & no.ori) {
      genes.with.tx <- unique(subset(net.df, Target == target | Ori_tx == target)$Gene)
    } else if (no.dup) {
      genes.with.tx <- unique(subset(net.df, Target == target)$Gene)
    } else {
      stop("; Specify at least one of the constraints: no.dup = TRUE or no.ori = TRUE")
    }
    
    
    ## Detect if the target is 'problematic'
    if (length(genes.with.tx) >= half.set.size) {
      problematic.flag <- 1
      no.ori           <- FALSE
      genes.with.tx    <- unique(subset(net.df, Target == target)$Gene)
      message("; Target ", target, " is present in more than half of the networks. The 'no.ori == T' constraint is not applicable")
    }
    
    
    ## Return a set of genes without the target
    ## These are the genes that can receive the new target
    available.genes <- genes.set[!genes.set %in% genes.with.tx]
    
    ## These are the genes that can transfer a target
    # if (problematic.flag) {
    #   available.genes <- genes.set
    # }
    exchange.subset <- subset(net.df, Gene %in% available.genes & !Target %in% tx.in.net & !Target %in% ori.tx.in.net)
    
    if (!nrow(exchange.subset)) {
      message("; There is no a gene network to exchange targets with")
      return(NULL)
    }
    
    if (sum(exchange.subset$Realoc)) {
      exchange.subset <- subset(exchange.subset, Realoc == 1)
    }
    
    ## Choose an entry to exchange
    # to.entry <- subset(net.df, ID %in% exchange.subset)
    
    return(exchange.subset)
    
    
    ## When the Gene-Target association does not exist  
  } else {
    message("; Gene ", gene, " network does not have target ", target)
    return(NULL)
  }
}


###########
## Debug ##
###########
# results.dir <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/Test1/Rand_net"
# net.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/examples/data/Weighted_net_example.txt"
# net.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/Test1/Rand_net/TargetScan/TargetScan_clean_network.tab"
# nb.rand.net <- 1
# nb.cores <- 1
# exclude.targets <- T
# suffix.net.name <- "testing_suffix"


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument (Output directory): results.dir ")
  
} else if (!exists("net.file")) {
  stop("Missing mandatory argument: net.file ")
  
}

rand.net.dir <- file.path(results.dir, "Random_networks")
dir.create(rand.net.dir, showWarnings = F, recursive = T)



####################################
## Generate network (xseq format) ##
####################################
## Convert a data.frame network (Gene, Target, Weight) in a list of lists.
## Each sub-list is named as the Gene and their content correspond to their targets
## and the value of association (Weight).
df.net.to.xseq <- function(network){
  
  ## Rename columns
  colnames(network) <- c("Gene", "Target_gene", "Assoc_score")
  
  ##
  net.list <- split(network, f = network$Gene)
  
  xseq.net <- lapply(net.list, function(l){
    
    scores <- l$Assoc_score
    names(scores) <- l$Target_gene
    
    l <- scores
    
  })
  
  return(xseq.net)
  
}


#######################
## Read network file ##
#######################
message("; Reading network file")
net <- fread(net.file)
# net <- net[sample(1:nrow(net), size = 100000),]
# net <- net[1:100000,]

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
message("; Sorting network by descending number of targets")
net <- net %>% 
      group_by(Gene) %>% 
      mutate(Nb_target = n()) %>% 
      arrange(desc(Nb_target)) %>% 
      select(Gene, Target, Weight)

net.list  <- split(net, f = net$Gene)

## Sort the networks by decreasing size
gen.net.size.order <- order(unlist(lapply(net.list, nrow)), decreasing = T)
net.list           <- net.list[gen.net.size.order]


## The universe of target genes (many of them are repeated because are target of many genes)
## We keep this number to maintain exactly the same number of each target in the random network
all.targets <- as.vector(net$Target)

registerDoParallel(nb.cores)
new.targets.list <- NULL
message("; Generating random networks")
new.targets.list <- foreach(i = 1:nb.rand.net, .combine = 'cbind') %dopar% {

  lapply(net.list, function(l){

    message("; Number of targets to reallocate: ", length(all.targets))
    
    ## Network size/Number of targets
    nb.entries <- nrow(l)
    
    ## Shuffle the names on each iteration
    all.targets <<- sample(all.targets)
    
    ## This flag indicates whether the network should be modified after reallocation
    realoc.flag <- 0
    
    ## Original targets in the current network 
    ori.targets <- as.vector(l$Target)
    
    ## Get a new set of non-duplicated target genes
    ## 1) Get a set of entries that were not in the original network
    if (exclude.targets) {
      

      
      ## The remaining targets that are not part of the current network
      remaining.unique.names <- length(unique(all.targets[!all.targets %in% ori.targets]))
      
      # message("; Nb entries: ", nb.entries)
      # message("; Nb remaining unique: ", remaining.unique.names)
      
      if (nb.entries > remaining.unique.names) {
        
        new.genes <- all.targets[1:nb.entries]
        
        
        if (sum(ori.targets %in% new.genes)) {
          realoc.flag <- 1
        }
        
        
      } else {
        new.genes   <- unique(all.targets[!all.targets %in% ori.targets])[1:nb.entries]
      }

    ## 2) Get a random set of entries, non-duplicated but the randomly selected
    ## targets may overlap with the original network
    } else {
      new.genes <- unique(all.targets)[1:nb.entries]
    }
    
    
    
    ## Update the targets' vector
    if (realoc.flag) {
      all.targets  <<- all.targets[-(1:nb.entries)]
    } else {
      new.genes.ind <- match(new.genes, all.targets) 
      all.targets  <<- all.targets[-new.genes.ind]
    }

    
    data.frame(Gene   = unique(l$Gene),
               Target = new.genes,
               Realoc = realoc.flag,
               Ori_tx = ori.targets)
  })
}


## Convert the lists into a dataframe
rand.net.df        <- purrr::map_dfr(new.targets.list, data.table)
rand.net.df$Weight <- sample(net$Weight)
# rand.net.df.cp <- rand.net.df
#rand.net.df <- rand.net.df.cp


if (exclude.targets) {
  
  ## These 'problematic' genes are those genes targeted by more than the half of the 
  ## genes in the network, therefore in these cases the exclusive target condition
  ## cannot be applied.
  nb.genes.th    <- floor(length(unique(rand.net.df$Gene))/2)
  problematic.tx <- as.vector(unlist(rand.net.df %>% 
                                       group_by(Target) %>% 
                                       tally() %>% 
                                       arrange(desc(n)) %>% 
                                       dplyr::filter(n >= nb.genes.th) %>% 
                                       select(Target)))
  
  
  ## Clean the reallocation column
  ## First: get the positions where the random target is already in the original set of targets
  rand.net.df <- rand.net.df %>% 
    arrange(Gene)
  genes.order <- unique(rand.net.df$Gene)
  rand.net.list <- split(rand.net.df, f = rand.net.df$Gene)
  rand.net.list <- rand.net.list[genes.order]
  Realoc.tx.in.ori <- sapply(rand.net.list, function(l){
    
    ## Original targets
    targets <- as.vector(l$Ori_tx)
    
    ## Boolean: random targets in original set of targets?
    l$Target %in% targets
    
  })
  Realoc.tx.in.ori <- as.vector(unlist(Realoc.tx.in.ori))
  
  ## Second: get the positions where the target gene is duplicated within the same gene network
  Realoc.tx.dup <- sapply(rand.net.list, function(l){
    
    ## IS the random target duplicated ?
    duplicated(l$Target)
    
  })
  Realoc.tx.dup <- as.vector(unlist(Realoc.tx.dup))
  
  ## If at least one of this conditions is true, then Realoc value is 1
  clean.realoc.vector <- ifelse(Realoc.tx.in.ori + Realoc.tx.dup >= 1, yes = 1, no = Realoc.tx.in.ori + Realoc.tx.dup)
  rand.net.df$Realoc  <- clean.realoc.vector
  
  ## Add and ID to each entry
  rand.net.df$ID          <- paste0("Entry_", 1:nrow(rand.net.df))
  sum(clean.realoc.vector)
  
  
  ##############################################################################
  ##############################################################################
  
  rand.net.realoc <- subset(rand.net.df, Realoc == 1)
  no.output <- sapply(1:nrow(rand.net.realoc), function(l){
    
    to.realoc <- sum(rand.net.df$Realoc == 1)
    message("; Targets to reallocate: ", to.realoc)
    
    ## Run when is necessary to reallocate
    if (to.realoc) {
      entry.from <- rand.net.realoc[l,]$ID     ## This is the target that will be sent (from)
      
      ## Check that the entry is not correctly reallocated
      if (subset(rand.net.df, ID == entry.from)$Realoc) {
        
        gg <- rand.net.realoc[l,]$Gene
        tt <- rand.net.realoc[l,]$Target
        
        entry.exists <- subset(rand.net.df, Gene == gg & Target == tt)
        
        if (nrow(entry.exists)) {
          
          ## Find the available positions
          available.entries <- find.new.gene.target(net.df = rand.net.df,
                                                    gene   = gg,
                                                    target = tt,
                                                    no.dup = T,
                                                    no.ori = T)
          
          if (nrow(available.entries)) {
            
            ## This is the target that will be received (to)
            entry.to <- sample(available.entries$ID, 1)
            
            ## Exchange targets
            ## tt == rand.net.df$Target[which(rand.net.df$ID == entry.from)]
            new.tt <- rand.net.df$Target[which(rand.net.df$ID == entry.to)]
            
            if (tt == new.tt) {
              stop("; Target and new target are the same")
            }
            
            rand.net.df$Target[which(rand.net.df$ID == entry.to)]   <<- tt
            rand.net.df$Target[which(rand.net.df$ID == entry.from)] <<- new.tt
            
            rand.net.df$Realoc[which(rand.net.df$ID == entry.to)]   <<- 0
            rand.net.df$Realoc[which(rand.net.df$ID == entry.from)] <<- 0
            
            ## In case there are no available entries to exchange targets
          } else {
            message("; Entry: ", entry.from, " has no available entries to exchange")
            NULL
          }
        }
        
        ## In case the entry is already reallocated
      } else {
        message("; This entry is already correctly reallocated")
        NULL
      }
    } else {
      message("; No need for target reallocation")
    }
  })
}
# new.targets.df <- lapply(lapply(new.targets.list, lapply, data.frame), purrr::map_dfr, data.table)


rand.net.df <- rand.net.df %>% 
                select(Gene, Target, Weight)

## Rename the iteration number in the output file, only when the --suffix option
## is indicated by the user
if (!is.null(suffix.net.name)) {
  it <- suffix.net.name
} else {
  it <- 1
}

#####################################
## Export networks as Rdata object ##
#####################################
rand.net.list <- df.net.to.xseq(rand.net.df)
rand.net.rds <- file.path(rand.net.dir, paste0("Random_network_", it, ".rds"))
message("; Exporting random network as Rdata object: ", rand.net.rds)
saveRDS(rand.net.list, file = rand.net.rds)

#################################
## Export network as text file ##
#################################
rand.net.file <- file.path(rand.net.dir, paste0("Random_network_", it, ".tab"))
message("; Exporting random network as text file: ", rand.net.file)
fwrite(rand.net.df, file = rand.net.file, sep = "\t", row.names = F, col.names = T)


## To check

# ## Export the networks as text files
# it <- 0
# lapply(rand.net.df, function(l){
#   
#   colnames(l) <- c("Gene", "Target", "Weight")
#   
#   
#   ## Rename the iteration number in the output file, only when the --suffix option
#   ## is indicated by the user
#   if (!is.null(suffix.net.name)) {
#     it <<- suffix.net.name
#   } else {
#     it <<- it + 1
#   }
#   
#   
#   #####################################
#   ## Export networks as Rdata object ##
#   #####################################
#   rand.net.list <- df.net.to.xseq(l)
#   rand.net.rds <- file.path(rand.net.dir, paste0("Random_network_", it, ".rds"))
#   message("; Exporting random network as Rdata object: ", rand.net.rds)
#   saveRDS(rand.net.list, file = rand.net.rds)
#   
#   #################################
#   ## Export network as text file ##
#   #################################
#   rand.net.file <- file.path(rand.net.dir, paste0("Random_network_", it, ".tab"))
#   message("; Exporting random network as text file: ", rand.net.file)
#   fwrite(l, file = rand.net.file, sep = "\t", row.names = F, col.names = T)
# })