## Packages
required <- c("tidyverse", "ggpubr", "ggbiplot", "devtools", "gplots", "UsefulFunctions")
for (lib in required)
{
  if (!require(lib, character.only = T))
  {
    install.packages(lib)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

# Comparisons and Colours
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"), c("VD1.CD27HI", "CD8.EMRA"), c("VD1.CD27HI", "CD8.Naive"), c("VD1.CD27HI", "VD2"),
                       c("VD1.CD27LO", "CD8.EMRA"), c("VD1.CD27LO", "CD8.Naive"), c("VD1.CD27LO", "VD2"),
                       c("CD8.EMRA", "CD8.Naive"), c("CD8.EMRA", "VD2"), c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999", "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00", "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")


load("Bulk/Final_edgeR.RData")
getwd()
setwd("./Bulk/Output/Cyto/")

tryCatch(expr = { library("limma")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("limma")}, 
         finally = library("limma"))

tryCatch(expr = { library("GSA")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("GSA")}, 
         finally = library("GSA"))

tryCatch(expr = { library("RCurl")}, 
         error = function(e) { 
           install.packages("RCurl")}, 
         finally = library("RCurl"))

# gmt_url <- "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"
# 
# #list all the files on the server
# filenames <- getURL(gmt_url)
# tc <- textConnection(filenames)
# contents <- readLines(tc)
# close(tc)
# 
# #get the gmt that has all the pathways and does not include terms 
# #inferred from electronic annotations(IEA)
# #start with gmt file that has pathways only
# rx <- gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
#               contents, perl = T)
# gmt_file <- unlist(regmatches(contents, rx))
# working_dir <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/"
# 
# dest_gmt_file <- file.path(working_dir, paste("Supplementary_Table3_", gmt_file, sep="") )
# 
# download.file(
#   paste(gmt_url,gmt_file,sep=""),
#   destfile = dest_gmt_file
# )


# Geneset enrichment
## GO_genesets
# load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/GO_genesets.rdata")
library(qusage)
All_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/msigdb.v6.2.entrez.gmt")
KEGG_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/c2.cp.kegg.v6.2.entrez.gmt")
GO_terms <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/c5.all.v6.2.entrez.gmt")


# Filter 
## Genesets that only appear in KEGG or GO databases (6103 genesets)
filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_terms)]

## Genesets for very small/very big sizes (reduces multiple comparison deficit) (4326 genesets)
geneset_sizes <- unlist(lapply(filter_gmt, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- filter_gmt[geneset_indices]

# Perform camera analysis
idx <- ids2indices(filtered_set, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])


# BiocManager::install("qusage")
library(qusage)
camera_results_a <- camera_results[rownames(camera_results) %in% names(filtered_set),]

# sig_camera_results <- droplevels(subset(camera_results, FDR <= 0.01))
# GO_GD <- rownames_to_column(sig_camera_results, var = "GO_Pathway")

genesets_filtered <- idx
data_for_gs_analysis <- v

camera_descr <- unlist(lapply(rownames(camera_results_a), 
                              function(x){unlist(strsplit(x,"\\%"))[1]}))
camera_Phenotype <- unlist(lapply(camera_results_a[, "Direction"], 
                                  function(x){if(x=="Up"){1}else{(-1)}}))

camera_genes <- c()
for(i in 1:length(rownames(camera_results_a))){
  current_geneset <- unlist( 
    genesets_filtered[which(names(genesets_filtered) %in% 
                              rownames(camera_results_a)[i])])
  current_genes <- c()
  for(j in 1:length(current_geneset)){
    if(j==length(current_geneset)){
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]],
                             sep = "")
    } else {
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]], ",", 
                             sep = "")
    }
  }
  camera_genes <- rbind(camera_genes, current_genes)
}
rownames(camera_genes) <- rownames(camera_results_a)

head(camera_results_a)

camera_results_generic_em <- data.frame(rownames(camera_results_a), camera_descr, 
                                        PValue = camera_results_a[, "PValue"],
                                        FDR = camera_results_a[, "FDR"],
                                        camera_Phenotype,
                                        camera_genes)

camera_results_file <- "camera_results_generic_em.txt"
write.table(camera_results_generic_em, file.path(camera_results_file), 
            col.name = T, sep = "\t", row.names = F, quote = F)
expression_file <- "expression_file.txt"

exp_fil <- as.data.frame(v$E)
write.table(exp_fil, file.path(expression_file),
            col.name = T, sep = "\t", row.names = F, quote = F)



#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library(RCy3)}, 
         error = function(e) { install_github("cytoscape/RCy3")}, finally = library(RCy3))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_threshold <- "0.05"
qvalue_threshold <- "0.001"

similarity_threshold <- "0.25"
similarity_metric <- "JACCARD"

# generic_gmt_file <- file.path(getwd(), gmt_file)
analysis_name <- "CD27LO_vs_CD27HI"
cur_model_name <- paste("camera", analysis_name, sep="_")
results_filename <- file.path(getwd(),  camera_results_file)






#######################################
#create EM - camera results
#######################################
current_network_name <- paste(cur_model_name, pvalue_threshold, qvalue_threshold, sep = "_")


results_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/camera_results_generic_em.txt"
expression_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/expression_file.txt"
generic_gmt_file <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/c5.all.v6.2.entrez.gmt"
em_command = paste('enrichmentmap build analysisType="generic"',
                   'gmtFile=', generic_gmt_file,
                   'pvalue=', pvalue_threshold,
                   'qvalue=', qvalue_threshold,
                   'similaritycutoff=', similarity_threshold,
                   'coefficients=', similarity_metric,
                   'enrichmentsDataset1=', results_filename,
                   'expressionDataset1=', expression_filename)

# Above had sep = " "

#enrichment map command will return the suid of newly created network.
# install.packages("BiocManager")
# BiocManager::install("RCy3")
library(RCy3)
response <- commandsGET(em_command)

current_network_suid <- 0
#enrichment map command will return the suid of newly created network unless it Failed.  
#If it failed it will contain the word failed
if(grepl(pattern="Failed", response)){
  paste(response)
} else {
  current_network_suid <- response
}




response <- renameNetwork(current_network_name, as.numeric(current_network_suid))
