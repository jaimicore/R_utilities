#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "GenomicRanges",
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
              help = "Output directory to export the tables (Mandatory)", metavar = "character"),
  
  make_option(c("-t", "--template"), type = "character", default = "NULL", 
              help = "File used as template to generate random positions. (Mandatory) ", metavar = "character"),  
  
  make_option(c("-f", "--format"), type = "character", default = "NULL", 
              help = "Input file format with genomic coordinates [BED | VCF] (Mandatory) ", metavar = "character"),  
  
  make_option(c("-s", "--side"), type = "character", default = "both", 
              help = "Insert the positions randomly in [left] side, [right] side or [both].", metavar = "character"),  
  
  make_option(c("-m", "--max_distance"), type = "numeric", default = 10, 
              help = "Max distance where the random mutations are generated, relative to the original coordinates in the template. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-n", "--min_distance"), type = "numeric", default = 10, 
              help = "Max distance where the random mutations are generated, relative to the original coordinates in the template. [Default \"%default\"] ", metavar = "number"),

  make_option(c("-i", "--new_ids"), type = "numeric", default = 0, 
              help = "Generate unique identifiers based on the given column. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-p", "--id_prefix"), type = "character", default = "Sample", 
              help = "Prefix used to rename the new IDs. This prefix will be followed by a number [Default \"%default\"] ", metavar = "character")
  
  );
message("; Reading arguments from command line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);



########################
## Set variable names ##
########################
results.dir         <- opt$output_directory
template.tab.file   <- opt$template
input.format        <- tolower(opt$format)
insert.side         <- opt$side
max.distance        <- as.numeric(opt$max_distance)
min.distance        <- as.numeric(opt$min_distance)
new.ids.col         <- as.numeric(opt$new_ids)
id.prefix           <- opt$id_prefix


#########################
## Mandatory variables ##
#########################
if (!exists("results.dir")) {
  stop("Missing mandatory argument (Output directory): results.dir ")
  
} else if (!exists("template.tab.file")) {
  stop("Missing mandatory argument: template.tab.file ")
  
} else if (!exists("input.format")) {
  stop("Missing mandatory argument: input.format ")
  
}
dir.create(results.dir, showWarnings = F, recursive = T)


###############################
## Generate random positions ##
###############################
## Position: numeric vector
## Max distance : max distance randomly generated from the provided position(s)
## Min distance : max distance randomly generated from the provided position(s)
## Side         : the new position should be at left, right or both sides [left|right|both].
## Iter         : number of times the function will try to randomly generate a non-negative integer
runif.positon <- function(position = 25,
                          max.distance = 10,
                          min.distance = 1,
                          side = "both",
                          iter = 5){
  
  na.pos <- 1
  max.iter <- 0
  while (na.pos == 1 & max.iter <= iter) {
    
    max.iter <- max.iter + 1
    
    ## Generate a new random position from a uniform distribution
    ## with limits in the max and min distance relative the position
    shift.pos <- round( runif(1, min = min.distance, max = max.distance) )
    
    
    ## Add direction to the shifted position
    ##
    ## Left : multiply by -1
    ## Right: multiply by +1 (left as it it)
    ## Both : randomly choose -1|+1 and multiply the shift
    if (side == "left") {
      
      shift.pos <- shift.pos * -1
      
    } else if (side == "right") {
      
      shift.pos <- shift.pos * 1
      
    } else if (side == "both") {
      
      shift.pos <- shift.pos * sample(c(-1,1), size = 1, replace = T)
    } else {
      stop("; 'side' parameter not recognized. Options: left|right|both")
    }
    
    
    ## Check that new position is not below to 0
    new.position <- position + shift.pos
    new.position <- ifelse(new.position >= 1, yes = new.position, no = NA)
    
    if (!is.na(new.position)) {
      na.pos <- 0
    }
  }
  
  
  ## Return a warning when after the N iterations the shifted position is lower than 1
  if (is.na(new.position)) {
    warning("; The shifted position of ", position, " is lower than 1 after ", iter," tries: NA returned.")
  }

  
  return(new.position)
}



#################################
## Read input coordinates file ##
#################################
message("; Reading template file: ", template.tab.file)
template.tab <- fread(template.tab.file)
template.header <- colnames(template.tab)

## Select the column from which new IDs will be generated
new.id.col.data <- template.tab[, ..new.ids.col]

## Separate the coordinates and supplementary columns
if (input.format == "vcf") {
  
  ## Col 1: chr
  ## Col 2: position
  template.coordinates <- template.tab[,1:2]
  colnames(template.coordinates) <- paste0("V",1:2)
  template.supp <- template.tab[,-c(1,2)]
  
} else  if (input.format == "bed") {
  
  ## Col 1: chr
  ## Col 2: start
  ## Col 3: end
  template.coordinates <- template.tab[,1:3]
  colnames(template.coordinates) <- paste0("V",1:3)
  template.supp <- template.tab[,-c(1,2,3)]
  
  gap.bed <- template.coordinates$V3 - template.coordinates$V2
  
}


###############################
## Generate random positions ##
###############################
message("; Generating random positions from a uniform distribution")
shifted.coordinates.start <- vapply(as.vector(unlist(template.coordinates[,2])),
                                    runif.positon, max.distance = 15, min.distance = 2, side = "both",
                                    numeric(1))

## Generate shifted coordinates
new.coordinates <- template.coordinates
new.coordinates$V2 <- shifted.coordinates.start

if (input.format == "bed") {
  
  if (length(new.coordinates$V2) == length(gap.bed) ) {
    new.coordinates$V3 <- new.coordinates$V2 + gap.bed
  } else {
    stop("; Error: the number of new coordinates and the length of the vector of gaps is not the same")
  }
  
}


## Concat new coordinates with supp columns
coords <- cbind(template.coordinates, template.supp)


########################################
## Generate random IDs (if indicated) ##
########################################
if (new.ids.col != 0) {
  
  old.ids <- unique(unlist(new.id.col.data))
  new.ids <- sample(paste0(id.prefix, sprintf("%04d", 1:length(old.ids))))
  
  ids.dict <- data.frame(Old = old.ids,
                         New = new.ids)
  
  ## Export dictionary
  IDs.dic.file <- file.path(results.dir, "Dictionary_IDs.txt")
  message("; Exporting ID dictionary: ", IDs.dic.file)
  fwrite(ids.dict, file = IDs.dic.file, sep = "\t", row.names = F, col.names = T)
  
  coords <- merge(x = coords, y = ids.dict, 
                  by.x = names(coords)[new.ids.col],
                  by.y = "Old")
  
  coords <- coords[,-1]
  
  if (input.format == "bed") {
    
    colnames(coords)[1:3] <- template.header[1:3]
    coords <- GenomicRanges::sort(GRanges(coords))
    coords <- data.frame(coords)[,-c(4,5)]
    colnames(coords) <- template.header
    
  } else if (input.format == "vcf") {
    
    colnames(coords) <- template.header
    
  }

}
colnames(coords) <- template.header


###############################################
## Remove positions with NAs and export them ##
###############################################
message("; Removing positions that were not shifted")
keep.positions <- which(complete.cases(coords))
remove.positions <- which(!complete.cases(coords))

removed.bed <- template.coordinates[remove.positions,]
coords <- coords[keep.positions,]

base.input <- basename(template.tab.file)
base.input <- unlist(strsplit(base.input, split = "\\."))

## Export coordinates file with positions that were not modified
nonshift.bed.file <- file.path(results.dir, paste0(base.input[1], "_non_shifted_positions.", base.input[2]))
message("; Exporting coordinates file with non-shifted positions: ", nonshift.bed.file)
fwrite(removed.bed, file = nonshift.bed.file, sep = "\t", row.names = F, col.names = F)


####################################################
## Export coordinates file with shifted positions ##
####################################################
coords.file <- file.path(results.dir, paste0(base.input[1], "_randomly_shifted.", base.input[2]))
message("; Exporting coordinates file with shifted positions: ", coords.file)
fwrite(coords, file = coords.file, sep = "\t", row.names = F, col.names = T)