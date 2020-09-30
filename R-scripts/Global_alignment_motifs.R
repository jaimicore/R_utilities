## Computes the DNAFULL substitution matrix
## considering all letters in the extended IUPAC consensus format
DNA_subs_matrix <- function() {
  
  ## Taken from: http://rosalind.info/glossary/dnafull/
  
  ## DNA consensus characters
  consensus.char <- c("A", "T", "G", "C", "S", "W", "R", "Y", "K", "M", "B", "V", "H", "D", "N")
  
  ##       A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
  A1 <- c( 5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2)  ## A
  T1 <- c(-4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2)  ## T
  G1 <- c(-4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2)  ## G
  C1 <- c(-4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2)  ## C
  S1 <- c(-4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1)  ## S: [C,G]   Strong interactions 
  W1 <- c( 1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1)  ## W: [A,T]   Weak interactions
  R1 <- c( 1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1 ,-1)  ## R: [A,G]   puRine
  Y1 <- c(-4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1)  ## Y: [C,T]   pYrimidine
  K1 <- c(-4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1)  ## K: [G,T]   Keto group
  M1 <- c( 1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1)  ## M: [A,C]   aMino group
  B1 <- c(-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1)  ## B: [C,G,T] Not A
  V1 <- c(-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1)  ## V: [A,C,G] Not T
  H1 <- c(-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1)  ## H: [A,C,T] Not G
  D1 <- c(-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1)  ## D: [A,G,T] Not C
  N1 <- c(-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)  ## N: aNy
  
  DNA.sub.mat <- rbind(A1, T1, G1, C1, S1, W1, R1, Y1, K1, M1, B1, V1, H1, D1, N1)
  colnames(DNA.sub.mat) <- consensus.char
  rownames(DNA.sub.mat) <- consensus.char

  return(DNA.sub.mat)
}



## Computes the score via backtracking the best scores in the matrix
backtrack_global <- function(seq1, seq2, backtrack, subs_matrix, int_penalty, ext_penalty){
  
  m <- length(seq1)
  n <- length(seq2)
  
  gap.counter <- 0
  score       <- 0
  align_seq1  <- NULL
  align_seq2  <- NULL
  bt.matrix   <- matrix('', m + 1, n + 1, dimnames = list(c('', seq1),c('', seq2)))
  
  ## Starts from the lower-right entry in the backtrack matrix
  
  while (m > 0 | n > 0) {
    
    bt.matrix[m + 1, n + 1] <- backtrack[m + 1, n + 1]
    
    #################
    ## Move up '|' ##
    #################
    if (backtrack[m + 1, n + 1] == '|') {
      
      gap.counter <- gap.counter + 1
      
      ## When moving down '|', add a gap in sequence 2
      align_seq1 <- c(align_seq1, seq1[m])
      align_seq2 <- c(align_seq2, '-')
      
      ## Penalty for opening a gap
      score <- score + (int_penalty * (gap.counter * 0.5))
      
      ## Jumps one row above 
      m <- m - 1
      
    
    ####################### 
    ## Move to left '--' ##
    #######################
    } else if (backtrack[m + 1, n + 1] == '--') {
      
      gap.counter <- gap.counter + 1
      
      ## When moving to right '--', add a gap in sequence 1
      align_seq1 <- c(align_seq1, '-')
      align_seq2 <- c(align_seq2, seq2[n])
      
      ## Penalty for opening a gap
      score <- score + (int_penalty * gap.counter)
      
      ## Jumps one column to the left
      n <- n - 1
      
    ###################  
    ## Diagonal '\\' ##
    ###################
    } else {
      
      ## Add the corresponding nucleotides
      ## No gaps added in diagonal
      align_seq1 <- c(align_seq1, seq1[m]) 
      align_seq2 <- c(align_seq2, seq2[n])
      
      ## Add score corresponding to a match
      score <- score + subs_matrix[seq1[m], seq2[n]]
      
      ## Jumps to the upper-left position
      n <- n - 1
      m <- m - 1
    }
  }
  
  alignment <- c(paste0(rev(align_seq1), collapse = ''), paste0(rev(align_seq2), collapse = ''))
  
  ## Calculates the number of external gaps
  al1 <- gsub(alignment[1], pattern = "([A-Za-z])-([A-Za-z])", replacement = "\\1N\\2")
  al2 <- gsub(alignment[2], pattern = "([A-Za-z])-([A-Za-z])", replacement = "\\1N\\2")
  al1 <- gsub(al1, pattern = "[A-Za-z]", replacement = "")
  al2 <- gsub(al2, pattern = "[A-Za-z]", replacement = "")
  
  ## Adjust the final score
  ext_gap <- nchar(paste0(al1, al2)) * (int_penalty - ext_penalty)
  score <- score - ext_gap
  
  return(list(score     = score,
              alignment = alignment,
              path      = bt.matrix))
}


## Min-max normalization
## Input: a numeric vector
mm_normalise <- function(x, na.rm = TRUE) {
  ranx <- range(x, na.rm = na.rm)
  (x - ranx[1]) / diff(ranx)
}


library(compiler)
global_alignment <- cmpfun(function(seq1 = "GTTGCCATGGCAAC",
                                    seq2 = "GTTGCCAGGCAAC",
                                    subs_matrix = NULL,
                                    int_penalty = -5,
                                    ext_penalty = -2) {
  
  ## String to vector
  seq1 <- unlist(strsplit(seq1, ''))
  seq2 <- unlist(strsplit(seq2, ''))
  
  m <- length(seq1)
  n <- length(seq2)
  
  ## |  : Move down
  ## -- : Move to the right
  ## \\ : Move diagonally
  backtrack_key <- c('|', '--', '\\')
  
  ## Initialize matrix full of 0s and with an extra row and columns
  scoring.matrix <- matrix(0, length(seq1) + 1, length(seq2) + 1, dimnames = list(c('', seq1), c('', seq2)))
  
  ## Initialize first row and column by adding the gap score on each position
  scoring.matrix[1,] <- cumsum(c(0, rep(int_penalty, ncol(scoring.matrix) - 1)))
  scoring.matrix[,1] <- cumsum(c(0, rep(int_penalty, nrow(scoring.matrix) - 1)))
  
  
  ## Initialize backtrack matrix (visualization purposes)
  backtrack.matrix <- matrix('', length(seq1) + 1, length(seq2) + 1, dimnames = list(c('', seq1),c('', seq2)))
  backtrack.matrix[1,] <- '--'
  backtrack.matrix[,1] <- '|'
  backtrack.matrix[1,1] <- '\\'
  
  gap.counter <- 1
  for (i in seq(2, m + 1)) {
    for (j in seq(2, n + 1)) {
      
      ## Calculate the three possibles scores
      s_up   <- scoring.matrix[i - 1, j]     + int_penalty * gap.counter
      s_left <- scoring.matrix[i ,j - 1]     + int_penalty * gap.counter
      s_diag <- scoring.matrix[i - 1, j - 1] + subs_matrix[seq1[i - 1], seq2[j - 1]]
      
      scores <- c(s_up, s_left, s_diag) 
      
      ## Get the max score and fill the scoring matrix
      backtrack_update <- which.max(scores)
      score_matrix_update <- max(scores)
      
      # if (backtrack_update <= 3) {
      #   gap.counter <- gap.counter + 1
      # }
      
      ## Update the backtrack matrix with the corresponding symbol
      scoring.matrix[i,j]   <- score_matrix_update
      backtrack.matrix[i,j] <- backtrack_key[backtrack_update]
      
    }
  }
  
  backtrack <- backtrack_global(seq1,
                                seq2,
                                backtrack.matrix,
                                subs_matrix,
                                int_penalty,
                                ext_penalty)
  
  ## The score es normalized by the range of possible score values 
  score.range <- range(scoring.matrix)
  score.norm  <- mm_normalise(c(score.range,backtrack$score))[3]
  
  ## Aggregate score
  score.agg <- backtrack$score * score.norm
  
  return(list(Score     = backtrack$score,
              N_score   = score.norm,
              Agg_score = score.agg,
              Alignment = backtrack$alignment,
              Matrix    = scoring.matrix,
              Backtrack = backtrack.matrix,
              Path      = backtrack$path))
  
})





score_matrix_heatmap <- function(score_matrix    = NULL,
                                 backtrack_matrix = NULL,
                                 logo1           = NULL,
                                 logo2           = NULL) {

  library(ggplot2)
  library(reshape2)
  
  seq1 <- colnames(score_matrix)
  seq2 <- rownames(score_matrix)
  
  ## Assign unique names to the columns, this is required for melt the dataframe
  colnames(score_matrix) <- paste0(seq1, seq_along(seq1))
  rownames(score_matrix) <- paste0(seq2, seq_along(seq2))
  
  score_matrix_m           <- reshape2::melt(score_matrix)
  colnames(score_matrix_m) <- c("Seq2", "Seq1", "Score")
  
  backtrack_matrix_m           <- reshape2::melt(backtrack_matrix)
  colnames(backtrack_matrix_m) <- c("Seq2", "Seq1", "Step")
  
  
  score_matrix_m <- cbind(score_matrix_m, Step = backtrack_matrix_m[,"Step"])
  
  ## Generate plot with scores
  lim.values <- max(abs(score_matrix))
  align.score.heatmap <- ggplot(score_matrix_m, aes(x = Seq1, y = rev(Seq2), fill = Score)) +
                          geom_tile() +
                          geom_text(aes(label = Step), colour = "white", size = 8) +
                          coord_equal() +
                          # scale_fill_distiller(palette = "RdPu", direction = +1, limits = c(-lim.values, lim.values)) +
                          # scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim.values, lim.values)) +
                          scale_fill_viridis_c(option = "A", direction = -1, limits = c(-lim.values, lim.values)) +
                          labs(x = "Sequence 1", y = "Sequence 2", title = "") +
                          scale_x_discrete(labels = seq1, position = "top") +
                          scale_y_discrete(labels = rev(seq2)) +
                          theme_classic()
  
  return(align.score.heatmap)
}


multiResultClass <- function(result1 = NULL, result2 = NULL)
{
  me <- list(
    Align_info = result1,
    Heatmap    = result2
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}



##########
## Main ##
##########
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

subs.mat <- DNA_subs_matrix()


TF.consensus.file         <- "/storage/mathelierarea/processed/jamondra/Projects/Global_motif_alignment/Motif_names.txt"
# TF.consensus.file         <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/Test1/Motif_names.txt"
TF.consensus              <- fread(TF.consensus.file, header = F)
colnames(TF.consensus)    <- c("TF", "Consensus", "Consensus_RC")
TF.consensus$Consensus    <- toupper(TF.consensus$Consensus)
TF.consensus$Consensus_RC <- toupper(TF.consensus$Consensus_RC)
TF.consensus$ID <- paste(TF.consensus$TF, 1:nrow(TF.consensus), sep = "_")
# View(TF.consensus)


TF.ref.name    <- "ZNF454_1126"
# TF.ref.name    <- "TFAP2A_759"
TF.ref.cons    <- subset(TF.consensus, ID == TF.ref.name)$Consensus
TF.ref.cons.rc <- subset(TF.consensus, ID == TF.ref.name)$Consensus_RC


comparison.tab  <- data.frame()
heatmaps.list   <- list()
int.penalty.gap <- -8
ext.penalty.gap <- -4
# for (tf in TF.consensus$ID) {

numCores <- detectCores() - 1

if (numCores > 50) {
  numCores <- 50
} else if (numCores <= 0) {
  numCores <- 1
}

registerDoParallel(numCores)
combined.output <- foreach(tf = TF.consensus$ID) %dopar% {
 
  result <- multiResultClass()
  message("; Comparing ", TF.ref.name, " vs ", tf)
  
  
  TF.target.cons    <- subset(TF.consensus, ID == tf)$Consensus
  TF.target.cons.rc <- subset(TF.consensus, ID == tf)$Consensus_RC
  
  ####################################
  ## Alignment in both orientations ##
  ####################################
  alig.info.F <- global_alignment(seq1        = TF.ref.cons,
                                  seq2        = TF.target.cons,
                                  subs_matrix = subs.mat,
                                  int_penalty = int.penalty.gap,
                                  ext_penalty = ext.penalty.gap)
  
  alig.info.R <- global_alignment(seq1        = TF.ref.cons,
                                  seq2        = TF.target.cons.rc,
                                  subs_matrix = subs.mat,
                                  int_penalty = int.penalty.gap,
                                  ext_penalty = ext.penalty.gap)
  
  
  #############################
  ## Choose best orientation ##
  #############################
  strand <- ""
  if (alig.info.F$N_score > alig.info.R$N_score ) {
    alig.info <- alig.info.F
    strand <- "F"
  } else if (alig.info.F$N_score < alig.info.R$N_score) {
    alig.info <- alig.info.R
    strand <- "R"
  } else {
    alig.info <- alig.info.F
    strand <- "F"
  }
  

  
  TF.compa.summary <- data.frame(Seq1_name  = TF.ref.name,
                                 Seq2_name  = tf,
                                 Seq1       = TF.ref.cons,
                                 Seq2       = TF.target.cons,
                                 Seq1_align = alig.info$Alignment[1],
                                 Seq2_align = alig.info$Alignment[2],
                                 Score      = alig.info$Score,
                                 NScore     = alig.info$N_score,
                                 Agg_Score  = alig.info$Agg_score,
                                 Strand     = strand)
  
  comparison.tab <- rbind(comparison.tab, TF.compa.summary)
  result$Align_info <- comparison.tab
  
  heatmaps.list[[TF.ref.name]][[tf]] <<- score_matrix_heatmap(score_matrix     = alig.info$Matrix,
                                                              backtrack_matrix = alig.info$Path)
  result$Heatmap <- heatmaps.list[[TF.ref.name]][[tf]]
  
  return(result)
  
}


## Create a DF with all the alignments
alignment.info <- lapply(combined.output, function(l){
                    l$Align_info
                  })
comparison.tab <- do.call(rbind, alignment.info)
comparison.tab <- comparison.tab %>% 
  arrange(desc(Agg_Score)) %>% 
  mutate(N_Agg_score = mm_normalise(Agg_Score))


## Create a list with all heatmaps
heatmaps.list <- list()
alignment.info <- lapply(combined.output, function(l){
                    m1 <- as.vector(l$Align_info$Seq1_name)
                    m2 <- as.vector(l$Align_info$Seq2_name)
                    # print(m1)
                    heatmaps.list[[m1]][[m2]] <<- l$Heatmap
                   })



ex.tf <- "SOX7_1699"
ex.tf <- "REST_647"
ex.tf <- "ZNF454_1127"
ex.tf <- "ATF4_130"
ex.tf <- "TFAP2A(var.2)_760"
ex.tf <- "TFAP2A(var.3)_761"
heatmaps.list[[TF.ref.name]][[ex.tf]]
View(comparison.tab[,-1])


##################################
## TF consensus table - Example ##
##################################

# convert-matrix -i Vertebrate_motifs_jaspar_2022.txt  -from  tf -to tf -return counts,consensus > Vertebrate_motifs_jaspar_2022_w_cons.txt
# cat Vertebrate_motifs_jaspar_2022_w_cons.txt | grep '^ID' | cut -d' ' -f3 > Motif_names.txt
# cat Vertebrate_motifs_jaspar_2022_w_cons.txt | grep 'IUPAC:'  | cut -d ' '  -f4 > Consensus_F.txt
# cat Vertebrate_motifs_jaspar_2022_w_cons.txt | grep 'IUPAC.rc:'  | cut -d ' '  -f4 > Consensus_R.txt



