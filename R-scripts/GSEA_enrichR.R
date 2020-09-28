#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "DT",
                        "enrichR",
                        "ggplot2",
                        "ggrepel",
                        "htmlwidgets",
                        "optparse",
                        "plotly")

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
  
  make_option(c("-l", "--gene_list"), type = "character", default = "NULL", 
              help = "A text file with a gene name per line. (Mandatory) ", metavar = "character"),  
  
  make_option(c("-n", "--N_most_enriched_terms"), type = "numeric", default = 20, 
              help = "Top N most enriched terms to display in the barplots. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-p", "--Pvalue_enrichment"), type = "numeric", default = 0.001, 
              help = "P-value to select relevant terms. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-g", "--geneset_name"), type = "character", default = "Cohort1", 
              help = "Cohort name. [Default \"%default\"] ", metavar = "character"),
  
  make_option(c("-t", "--title"), type = "character", default = "enrichR_analysis", 
              help = "Suffix for output files. [Default \"%default\"] ", metavar = "character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


########################
## Set variable names ##
########################
results.dir    <- opt$output_directory
gene.list.file <- opt$gene_list
geneset.name   <- opt$geneset_name
top.n          <- as.numeric(opt$N_most_enriched_terms)
enrich.pval    <- as.numeric(opt$Pvalue_enrichment)


## Example
# top.n          <- 20
# gene.list.file <- "../examples/data/dysregulated_genes.txt"
# results.dir    <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/R_utilities/examples/results/enrichR_results"
# enrich.pval    <- 0.001
# geneset.name   <- "Example_enrichR"


#########################
## Mandatory variables ##
#########################
if (!exists("gene.list.file")) {
  stop("Missing mandatory argument: gene.list.file ")
  
} else if (!exists("results.dir")) {
  stop("Missing mandatory argument: results.dir ")
  
}


###########################
## Create results folder ##
###########################
message("; Creating result folders")
out.folders <- list()
out.folders[["plots"]] <- file.path(results.dir, "plots")
out.folders[["tables"]] <- file.path(results.dir, "tables")
# out.folders[["RData"]] <- file.path(results.dir, "RData")
no.message <- sapply(out.folders, dir.create, recursive = TRUE, showWarnings = FALSE)


#################################
## Load filtered network table ##
#################################
message("; Reading gene list file: ", gene.list.file)
gene.list <- fread(gene.list.file, sep = "\t", header = FALSE)
colnames(gene.list) <- "Gene_names"
gene.list.vector <- unique(as.vector(gene.list$Gene_names))


##################
## Call enrichR ##
##################

## Databases
dbs <- c("KEGG_2019_Human",
         "WikiPathways_2019_Human",
         "Panther_2016",
         "GO_Biological_Process_2018")

## Run enrichR function
enrichr.tab.list <- enrichR::enrichr(genes = gene.list.vector,
                                     databases = dbs)

## Convert the list into a dataframe
enrichr.tab <- do.call(rbind, enrichr.tab.list)

## Add DB column
enrichr.tab$DB <- gsub(rownames(enrichr.tab), pattern = "\\.\\d+$", replacement = "", perl = T)

## Add Significance and cohort columns
enrichr.tab <- 
  enrichr.tab %>% 
  dplyr::filter(!grepl(Term, pattern = "Mus musculus")) %>% 
  mutate(Significance = round(-log10(Adjusted.P.value), digits = 2),
         Cohort = geneset.name,
         Nb_genes = lengths(gregexpr(';', Genes)) + 1 ) %>% 
  group_by(DB) %>% 
  arrange(desc(Significance), .by_group = T)

## Remove the GO ID
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = "\\(GO:\\d+\\)", replacement = "")
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = "Homo sapiens", replacement = "")
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = "__", replacement = "_")
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = " WP\\d+$", replacement = "", perl = T)
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = "_P0s\\d+$", replacement = "", perl = T)
enrichr.tab$Term <- gsub(enrichr.tab$Term, pattern = "P\\d+$", replacement = "", perl = T)


################################
## Ranking of enriched terms  ##
################################
GSEA.enrich.rank.plots <- list()
for (db in unique(enrichr.tab$DB)) {
  
  ## Filter the table by database
  enrichr.tab.db <- enrichr.tab %>% 
    dplyr::filter(DB == db) %>% 
    na.omit() %>% 
    arrange(desc(Significance))  
  
  enrichr.tab.db$Ranks <- 1:nrow(enrichr.tab.db)
  enrichr.tab.db$Lab <- enrichr.tab.db$Term
  enrichr.tab.db$Lab[top.n:nrow(enrichr.tab.db)] <- ""
  
  ## Generate plot
  enrichment.rank.plot <- NULL
  enrichment.rank.plot <- 
    ggplot(enrichr.tab.db, aes(x = Ranks, y = Significance, color = Nb_genes, label = Term)) +
    geom_line(data = enrichr.tab.db, aes(x = Ranks, y = Significance), inherit.aes = F, color = '#9ecae1', alpha = 0.45) +
    geom_point() +
    theme_bw() +
    xlab("Ranks") +
    ylab("Significance : -log10(p-value)") +
    geom_hline(aes(yintercept = -log10(enrich.pval)), linetype = "dashed", alpha = 0.5) +
    scale_colour_gradient(low = "#d0d1e6", high = "#023858") +
    ggtitle(paste0(db, " terms enrichment in ", geneset.name))
  
  GSEA.enrich.rank.plots[[db]][[geneset.name]] <- enrichment.rank.plot
  
  #############################
  ## Export interactive plot ##
  #############################
  enrichment.rank.plotly <- ggplotly(enrichment.rank.plot,
                                     tooltip = c("Ranks", "Significance", "Nb_genes", "Term")) %>%
                              config(displaylogo = FALSE,
                                     displayModeBar = FALSE)
  htmlwidgets::saveWidget(widget = enrichment.rank.plotly,
                          file = file.path(out.folders[["plots"]], paste0("Rankplot_terms_", geneset.name, "_", db, ".html")))

  
  ## Add labels
  enrichment.rank.plot <- enrichment.rank.plot +
    geom_label_repel(aes(label = Lab),
                     box.padding = 0.1,
                     # point.padding = 0.55,
                     # direction = "both",
                     force = 15,
                     size = 3.5,
                     color = "black",
                     max.iter = 20000,
                     # nudge_x = 50,
                     # nudge_y = 0.75,
                     segment.color = '#d9d9d9')
  
  message("; Exporting rank plots with all terms in ", db, " - Gene list: ", geneset.name)
  ggsave(plot = enrichment.rank.plot,
         filename = file.path(out.folders[["plots"]], paste0("Rankplot_terms_", geneset.name, "_", db, ".pdf")),
         device = "pdf",
         width = 13, height = 8)
}
# GSEA.enrich.rank.plots.rdata <- file.path(out.folders[["RData"]], paste0("Ranking_enrichment_terms_plots_", geneset.name, ".RData"))
# save(GSEA.enrich.rank.plots, file = GSEA.enrich.rank.plots.rdata)


## Filter significant terms
enrichr.tab.filt <- enrichr.tab %>% 
  dplyr::filter(Significance >= -log10(enrich.pval))

## Draw barplots with enriched terms
GSEA.enrich.plots <- list()
if (nrow(enrichr.tab.filt) > 1) {
  
  for (db in unique(enrichr.tab.filt$DB) ) {
    
    ## Filter the table by database
    enrichr.tab.filt.db <- enrichr.tab.filt %>% 
      dplyr::filter(DB == db) %>% 
      top_n(top.n) %>% 
      na.omit()
    
    ## Check the number of enriched terms
    if (nrow(enrichr.tab.filt) > 1) {
      
      ###############
      ## Draw plot ##
      ###############
      message("; Generating barplots with enriched terms in ", db)
      enrichment.plot <- ggplot(enrichr.tab.filt.db) +
        geom_bar(aes(x = Term, y = Significance, fill = Nb_genes), stat = "identity") +
        theme_classic() +
        theme(axis.text.y = element_text(hjust = 1, size = 10)) +
        scale_x_discrete(limits = rev(as.vector(enrichr.tab.filt.db$Term))) +
        labs(title = paste0("Enriched terms: ", geneset.name, " - ", db), x = "", y = "Significance -log10(p-value)") +
        coord_flip() +
        scale_fill_distiller(palette = "Blues", direction = +1, limits = c(0, max(enrichr.tab.filt.db$Nb_genes)))
      
      ##################
      ## Export plots ##
      ##################
      message("; Exporting barplots with enriched terms in ", db, " - Gene list: ", geneset.name)
      GSEA.enrich.plots[[db]][[geneset.name]] <- enrichment.plot
      ggsave(plot = enrichment.plot,
             filename = file.path(out.folders[["plots"]], paste0("Enriched_terms_", geneset.name, "_", db, ".pdf")),
             device = "pdf",
             width = 13, height = 8)
      
    } else {
      message("; No significant terms for ", db)
    }
    
  }
  
} else {
  message("; No significant terms enriched in the databases")
}
# GSEA.enrich.plots.rdata <- file.path(out.folders[["RData"]], paste0("Enriched_terms_plots_", geneset.name, ".RData"))
# save(GSEA.enrich.plots, file = GSEA.enrich.plots.rdata)

## Export filtered table
enrichr.tab.list.filter.file <- file.path(out.folders[["tables"]], paste0("Enriched_terms_", geneset.name, "_", geneset.name, ".tab"))
write.table(enrichr.tab.filt, file = enrichr.tab.list.filter.file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

enrichr.tab.list.file <- file.path(out.folders[["tables"]], paste0("All_terms_", geneset.name, "_", geneset.name, ".tab"))
write.table(enrichr.tab, file = enrichr.tab.list.file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# enrichr.tab.list.rdata <- file.path(out.folders[["RData"]], paste0("All_terms_", geneset.name, "_", geneset.name, ".Rdata"))
# save(enrichr.tab, file = enrichr.tab.list.rdata)