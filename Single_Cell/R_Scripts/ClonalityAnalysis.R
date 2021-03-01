# A script to investigate the clonality of the cells gained from Mike Stubbington analyses
## Original Author: Xi Chen (pre Seurat Version 2.3.4)
## Co-Author: Jack McMurray (Updated to Seurat 2.3.4, and then to version 4.0.0)

### Clonality values for each cell were calculated by looking for the TraCeR-derived TCR sequences in the iRepertoire data (represented as percentages).
### It's not entirely obvious which metric we should use to define a single 'clonality score' for a cell given that there are two loci but they are not always detected by TraCeR in every cell. 
### In reality it probably doesn't make a huge amount of difference.
set.seed(2021) # Year of reanalysis

# Load up packages and make common variables
cbcols <- c("VD1.CD27LO" = "#999999",
            "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00",
            "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")

## Loading tools and data ##
#### We're using the [Seurat](http://satijalab.org/seurat/) packages for single-cell analysis
library(Seurat)
library(tidyverse)
setwd("./Single_Cell/")
load("RData/ClonalityAnalyses.RData")
# load("RData/SeuratAnalyses.RData")
# 
# 
# 
# 
# clonality_values <- read.table("./Data/cell_clonality_from_iRep.csv", sep = ",", row.names = 1, header = T)
# 
# ### 10 in this table indicates that a sequence wasn't detected by TraCer.
# ### A value of 0 shows that the sequence deteced by TraCeR wasn't in the iRepertoire data.
# 
# name_map <- read.csv("Data/name_map.txt", sep = "\t", row.names = 1)
# rownames(name_map) = paste('X', rownames(name_map), sep = "")
# name_map$sort_class <- ifelse((name_map$sort_class == "naive"), "VD1.CD27HI", "VD1.CD27LO")
# name_map <- name_map[rownames(name_map) %in% rownames(gd.data[[]]), ]
# name_map <- setNames(as.character(name_map$sort_class), rownames(name_map))
# gd.data <- AddMetaData(object = gd.data, metadata = name_map, col.name = "sort_class")
# 
# 
# gd.data@meta.data <- merge(gd.data@meta.data, clonality_values, by = 0) %>% column_to_rownames(., var = "Row.names")
# 
# FeaturePlot(gd.data, c("D_clonality", "G_clonality", "max_clonality", "mean_clonality"), reduction = "pca")
# 
# ### So this demonstrates that high clonality (as measured by iRepertoire) is only found in the EMRA-like population.
# ### For any comparisons involving expanded and non-expanded cells, we must only consider cells that are in the EMRA population to start with.
# ### Otherwise, any differences will be swamped by the EMRA vs naive differences which aren't interesting here.
# ### Can we divide cells into bins depending on their clonality values and then compare markers between these bins?
# 
# 
# 
# EMRA_clonality <- gd.data@meta.data[gd.data@meta.data$sort_class == "VD1.CD27LO", ]
# 
# head(EMRA_clonality)
# hist(EMRA_clonality$D_clonality)
# hist(EMRA_clonality$G_clonality)
# hist(EMRA_clonality$max_clonality)
# hist(EMRA_clonality$mean_clonality)
# 
# ### Let's divide the cells into two classes using max_clonality
# ### Low clonality (0-10%)
# ### High clonality (>20%)
# ### and see if we can identify discriminating markers between these populations.
# ### Ideally, this would be relatively robust to changes in these boundaries.
# ### We'll make a function to classify the cells and use this to annotate them for markers etc.
# ### The function will be applied to every row in the metadata table.
# ### First show if there are any cells at all in the naive sort class that were detected in the iRepertoire data.
# 
# classify_clonality <- function(x, low_ceiling = 10, high_floor = 20){
#     max_clonality = as.numeric(x["max_clonality"])
#     sort_class = x["sort_class"]
#     clone_class = "no_class"
#     #if (max_clonality>0 && sort_class == "EMRA") {
#     if (max_clonality>0) {
#         if (max_clonality < low_ceiling){
#             clone_class = "low_clonality"
#         } else if (max_clonality>high_floor){
#             clone_class = "high_clonality"
#         } else{
#         }
#     }
#     return(clone_class)
# }
# 
# apply(gd.data@meta.data, 1,
#       classify_clonality,
#       low_ceiling = 20,
#       high_floor = 20)
# 
# gd.data <- AddMetaData(gd.data, apply(gd.data@meta.data, 1, classify_clonality,
#                                       low_ceiling = 20, high_floor = 20), "clone_class_all")
# 
# ### Then do again for just the EMRA cells
# classify_clonality <- function(x, low_ceiling = 10, high_floor = 20){
#     max_clonality = as.numeric(x["max_clonality"])
#     sort_class = x["sort_class"]
#     clone_class = "no_class"
#     if (max_clonality>0 && sort_class == "VD1.CD27LO") {
#         if (max_clonality < low_ceiling){
#             clone_class = "low_clonality"
#         } else if (max_clonality>high_floor){
#             clone_class = "high_clonality"
#         } else{
#         }
#     }
#     return(clone_class)
# }
# 
# gd.data <- AddMetaData(gd.data, apply(gd.data@meta.data, 1, classify_clonality,
#                                       low_ceiling = 20, high_floor = 20),
#                        "clone_class_EMRA")
# 
# ### Below is the plot where we assign clonality classes to all the cells no matter their sort identity
# PCAPlot(gd.data, group.by = "clone_class_all")
# 
# ### Then we can compare this with the same plot if we only assign clonality classes to the cells sorted CD27lo cells.
# PCAPlot(gd.data, group.by = "clone_class_EMRA")
# 
# 
# ### So there are a couple in the naive sort that were at least detected in the iRepertoire data but really not many.
# ### Now re-do the EMRA clone class so that it contains the two bins discussed above.
# gd.data <- AddMetaData(gd.data, apply(gd.data@meta.data, 1, classify_clonality,
#                                       low_ceiling = 10, high_floor = 20), "clone_class_EMRA")
# PCAPlot(gd.data, group.by = "clone_class_EMRA")
# 
# 
# ### UNSURE?!
# gd.data_Clonality <- gd.data
# Idents(gd.data_Clonality) <- "clone_class_EMRA"
# 
# gd.data_Clonality@assays$RNA@scale.data <- as.matrix(gd.data_Clonality@assays$RNA@scale.data)
# 
# high_clonality_markers <- FindMarkers(gd.data_Clonality, ident.1 = "high_clonality", ident.2 = "no_class", min.pct = 0.25, only.pos = T)
# low_clonality_markers = FindMarkers(gd.data_Clonality, "low_clonality", "high_clonality", min.pct = 0.25, only.pos = T)
# 
# write.table(high_clonality_markers, file = "./Output/high_clonality_markers.txt", sep = "\t", quote = F)
# write.table(low_clonality_markers, file = "./Output/low_clonality_markers.txt", sep = "\t", quote = F)
# 
# ### Is there anything interesting going on with the set of markers plotted above when we look at high vs low clonality cells?
# VlnPlot(gd.data_Clonality,
#         c("PRF1",
#           "GZMA",
#           "CX3CR1",
#           "SELL",
#           "TCF7",
#           "LEF1",
#           "IL7R",
#           "LTB",
#           "GZMB"), ncol = 2)
# 
# high_clonality_markers[c("PRF1", "GZMA", "GZMB", "CX3CR1", "TCF7"), ]
# 
# ### Granzyme B looks to be significantly up in the high-clonality cells but perforin and GZMA not so much. Also, TCF7 appears higher in the high-clonality cells - that's a bit strange because that's typically associated with the naive cells.
# ### It's possible that these differences are actually due to a single large clone being different to the others. We should now look at differences between individual clones.
# #### I don't agree with this above... Seems to show the perfect story...
# clonotype_metadata <- read.table("./Data/clonotype_metadata.txt", row.names = 1)
# colnames(clonotype_metadata) = "clonotype"
# 
# gd.data_Clonality <- AddMetaData(gd.data_Clonality, clonotype_metadata[rownames(gd.data_Clonality@meta.data),,drop = F], "clonotype")
# head(gd.data_Clonality@meta.data)
# 
# ### A clonotype of 'X' means that it was a unique cell.
# #### Need to get a palette of 14 colours
# 
# # temp <- c(5,7,6,4,8)
# # barplot(temp, col = "#B53FE3") # test colours out
# 
# clonalityColours <- c("A" = "#DB724F", "B" = "#B53FE3", "C" = "#DAD65F", "D" = "#DD59B8",
#                       "E" = "#D7CC97", "F" = "#78D9DC", "G" = "#D2A0DF", "H" = "#766ED9",
#                       "I" = "#6FE18E", "J" = "#A0E2BC", "K" = "#94EA4E", "L" = "#79A0C9",
#                       "M" = "#D38B99", "X" = "#999999")
# 
# 
# cols_for_sort <- c("VD1.CD27LO" = "#999999", "VD1.CD27HI" = "#009E73")
# 
# pdf("./Figures/CellReportsRewrite/1F_Clonality.pdf")
# PCAPlot(gd.data_Clonality, group.by = "clonotype", cols = clonalityColours) +
#   theme(legend.title = element_blank(), legend.position = "top")
# dev.off()

# Is it statistically significant that EMRA-like cells have more cells in an expanded clone than not?
metaClonal <- gd.data_Clonality@meta.data
library(reshape2)

# Want to test if X matches with CD27HI exlcusively
metaClonal$expand <- ifelse((metaClonal$clonotype == "X"), "singlet", "expanded")
sorted <- metaClonal[, c("sort_class", "expand")]
continSort <- dcast(sorted, sort_class ~ expand) %>% column_to_rownames(., var  = "sort_class")

chisq_S <- chisq.test(continSort)

clustered <- metaClonal[, c("ident", "expand")]
continClust <- dcast(clustered, ident ~ expand) %>% column_to_rownames(., var  = "ident")
chisq_C <- chisq.test(continClust)




### So, are there markers that distinguish one clonotype from all the others? 
### Probably only worth doing for the larger clonotypes (A, E, J, K)

table(gd.data_Clonality@meta.data["clonotype"])
### So, clonotypes A, E, J and K worth looking at further I think


Idents(gd.data_Clonality) <- "clonotype"

## For each donor, find the largest clone
metaClonal %>% dcast(orig.ident ~ clonotype)

### Calculate markers for each individual clonotype
clonotype_A_markers <- FindMarkers(gd.data_Clonality, "A", c("E","J","K")) %>% rownames_to_column(., var = "SYMBOL")
clonotype_E_markers <- FindMarkers(gd.data_Clonality, "E", c("A","J","K")) %>% rownames_to_column(., var = "SYMBOL")
clonotype_J_markers <- FindMarkers(gd.data_Clonality, "J", c("A","E","K")) %>% rownames_to_column(., var = "SYMBOL")
clonotype_K_markers <- FindMarkers(gd.data_Clonality, "K", c("A","E","J")) %>% rownames_to_column(., var = "SYMBOL")

clonotype_JK_markers <- FindMarkers(gd.data_Clonality, c("J", "K"), c("A","E")) %>% rownames_to_column(., var = "SYMBOL")


# clonotype_A_markers <- FindMarkers(gd.data_Clonality, "A", c("E","K")) %>% rownames_to_column(., var = "SYMBOL")
# clonotype_E_markers <- FindMarkers(gd.data_Clonality, "E", c("A", "K")) %>% rownames_to_column(., var = "SYMBOL")
# clonotype_K_markers <- FindMarkers(gd.data_Clonality, "K", c("A", "E")) %>% rownames_to_column(., var = "SYMBOL")



write.csv(droplevels(subset(clonotype_A_markers, p_val_adj < 0.05)), file = "./Output/Clono_New/A_markers.csv", row.names = F)
write.csv(droplevels(subset(clonotype_E_markers, p_val_adj < 0.05)), file = "./Output/Clono_New/E_markers.csv", row.names = F)
write.csv(droplevels(subset(clonotype_K_markers, p_val_adj < 0.05)), file = "./Output/Clono_New/K_markers.csv", row.names = F)
write.csv(droplevels(subset(clonotype_J_markers, p_val_adj < 0.05)), file = "./Output/Clono_New/J_markers.csv", row.names = F)
write.csv(droplevels(subset(clonotype_JK_markers, p_val_adj < 0.05)), file = "./Output/Clono_New/JK_markers.csv", row.names = F)



### Since we only have four of interest we can also do the 2 vs 2 comparisons. Not sure how useful these will be.
### Since each group here contains half of the clonotype classes we can just do three comparisons and include both negative and positive markers. So, for clonotype_AE_markers the positive markers will be the ones up in both A and E whilst the negative markers will be the ones up in J and K.
# 
# clonotype_AE_markers = FindMarkers(gd.data_Clonality, c("A","E"), c("J","K"))
# clonotype_AJ_markers = FindMarkers(gd.data_Clonality, c("A","J"), c("E","K"))
# clonotype_AK_markers = FindMarkers(gd.data_Clonality, c("A","K"), c("E","J"))
# 
# ### Write out all the marker tables we've made here
# write.table(clonotype_A_markers, "Output/clonotype_A_markers.txt", sep = "\t", quote = F)
# write.table(clonotype_E_markers, "Output/clonotype_E_markers.txt", sep = "\t", quote = F)
# write.table(clonotype_J_markers, "Output/clonotype_J_markers.txt", sep = "\t", quote = F)
# write.table(clonotype_K_markers, "Output/clonotype_K_markers.txt", sep = "\t", quote = F)
# write.table(clonotype_AE_markers, "Output/clonotype_AE_markers.txt", sep = "\t", quote = F)
# write.table(clonotype_AJ_markers, "Output/clonotype_AJmarkers.txt", sep = "\t", quote = F)
# write.table(clonotype_AK_markers, "Output/clonotype_AK_markers.txt", sep = "\t", quote = F)


# Email from Carrie:
genes_for_violin = list(NK_receptors = c("KLRB1",
                                       "KLRG1",
                                       "KLRC2", "KLRC3", "KLRC4",
                                       "CADM1",
                                       "FCGR3A", "FCGR3B",
                                       "KIR2DL1", "KIR2DL2", "KIR2DL3",
                                       "KIR2DS2", "KIR3DS1",
                                       "KIR3DL2",
                                       "CD160",
                                       "KLRF1", "KLRK1", "KLRD1", "KLRC1",
                                       "LILRB1", "LILRB2"), 
                        cytokines_chemokines_and_receptors = c("CCL28",
                                                               "CCL3", "CCL4", "CCL5",
                                                               "TNFSF13B",
                                                               "TGFBR1",
                                                               "IL32",
                                                               "CXCR4",
                                                               "IL7R",
                                                               "CX3CR1"),
                        homing_adhesion = c("ITGA5", "ITGAD", "ITGAX", "ITGAE", "ITGAL",
                                            "SELL",
                                            "PECAM1",
                                            "S1PR5",
                                            "CD2",
                                            "SEMA4C", "SEMA4F", "SEMA3C", 
                                            "CD52",
                                            "CLEC11A"), 
                        effector_molecules = c("PRF1",
                                               "GZMA", "GZMB", "GZMK", "GZMM", "GZMH"),
                        transcription = c("NANOG", "TCF7", "PRDM1", "EOMES", "SOCS2"),
                        carrie = c("TRBV11-3", "TRGV3", "TRGV2", "TRGV5", "IL32", "GNLY", "KLRB1",
                                   "XIST", "FCGR3A", "CX3CR1"))

genes_for_violin

# pdf(file = "Output/test.pdf")
VlnPlot(gd.data_Clonality,c("KLRB1", "KLRG1", "KLRC2", "KLRC3"), group.by = "clonotype")
VlnPlot(gd.data_Clonality,c("KLRC4", "CADM1", "FCGR3A", "FCGR3B"), group.by = "clonotype")
dev.off()


save.image("./RData/ClonalityAnalyses.RData")
# load("final_clonality_analysis.RData")

vp = function(g, data_to_use = gd.data_Clonality){
  VlnPlot(data_to_use, g,
          group.by = "clonotype")
}

for (n in names(genes_for_violin)){ 
  print(n)
  file = paste("./Figures/CellReportsRewrite/Test/", n, ".pdf", sep="")
  pdf(file=file)
  genes = genes_for_violin[[n]]
  genes = split(genes, ceiling(seq_along(genes)/4))
  lapply(genes, vp)
  dev.off()
}

head(MD)

identifiers <- as.data.frame(gd.data_Clonality@ident)
head(identifiers)

carrie = c("TRGV3", "TRGV2", "TRGV5", "TRBV11-3",
           "XIST")

pdf("./Figures/CellReportsRewrite/Supp1B_Clone_diffs.pdf", height = 6, width = 10)
VlnPlot(gd.data_Clonality, carrie, 
        ncol = 5, idents = c("A", "E", "K", "J"), 
        cols = c("A" = "#000000", "E" = "#CC79A7", "K" = "#56B4E9", "J" = "#FF6600"))
dev.off()

## UNEEDED

# these_genes <- c("XIST", "TRGV3", "TRGV2", "TRGV5", "TRBV11-3")
# 
# tpms <- as.data.frame(gd.data_Clonality@data)
# tpm1 <-tpms %>% rownames_to_column(., var = "gene") %>% gather(key = "Cell", value = "tpm", -gene)
# tpm2 <- tpm1[tpm1$gene %in% these_genes,]
# clusts <- as.data.frame(gd.data_Clonality@ident) %>% rownames_to_column(., var = "Cell")
# colnames(clusts) <- c("Cell", "Cluster")
# 
# tpm3 <- merge(tpm2, clusts, by = "Cell")
# toplot <- tpm3[tpm3$Cluster %in% c("A", "E", "J"), ]
# 
# toplot$log_tpm <- log(toplot$tpm + 1)
# library(ggpubr)
# 
# dodge <- position_dodge(width = 0.9)
# 
# toplot$gene <- as.factor(toplot$gene)
# for(i in levels(toplot$gene)){
#   print(i)
#   working <- droplevels(subset(toplot, gene == i))
#   temp_plot <- ggplot(working, aes(x = Cluster, y = tpm, fill = Cluster)) +
#   stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.5) +
#   stat_summary(fun.y = mean, geom = "bar", position = dodge) + 
#   theme_bw() + 
#   # scale_fill_manual(values = cbcols) +
#   labs(x = "Clonotype", y = paste0("TPM_", i)) + 
#   stat_compare_means(label = "p.signif", comparisons = list((c("A", "E")), c("A", "J"), c("E", "J"))) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 12))
#   filen <- paste0(i, ".pdf")
#   ggsave(filen, plot = temp_plot, path = "./Single_Cell/Figures/Gamma-Delta/Supplementary")
#   }
# 
# getwd()
