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
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"), c("VD1.CD27HI", "CD8.EMRA"), c("VD1.CD27HI", "CD8.Naive"),
                       c("VD1.CD27HI", "VD2"), c("VD1.CD27LO", "CD8.EMRA"), c("VD1.CD27LO", "CD8.Naive"),
                       c("VD1.CD27LO", "VD2"), c("CD8.EMRA", "CD8.Naive"), c("CD8.EMRA", "VD2"),
                       c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999", "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00", "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")

# Soruce script
source("./Bulk/R_Scripts/Counts.R")

# Create copy of "x"
x1 <- x

# Figure 2A - MDS plot of all 
## Remove VD2
x1$samples <- droplevels(subset(x1$samples, group != "VD2"))
lcpm_noVD2 <- lcpm %>% as.data.frame() %>% dplyr:: select(-contains("VD2")) %>% as.matrix()

## Make colour scheme and plot
col.cell <- c("#56B4E9","#E69F00","#009E73","#999999")[x1$samples$group]
data.frame(x1$samples$group, col.cell)

pdf("../Figures/Figure_2/Figure2A.pdf", height = 6, width = 6)
plotMDS(lcpm_noVD2, pch = 16, cex = 2, col = col.cell, top = 14191)
legend("top",
       fill = c("#56B4E9",
                "#E69F00", "#009E73", "#999999"),
       legend = levels(x1$samples$group))
dev.off()


# Figure 2B - Volcano of differentially expressed genes in VD1.CD27HI vs VD1.CD27LO (annotate LTB, LEF1, and TCF7)
library(ggrepel)

## Get data in correct format
myData <- as.data.frame(CD27LO.vs.CD27HI)
myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)
myData$Labels <- myData$SYMBOL

labelled_genes <- c("LTB", "LEF1", "TCF7", "CCR7", "SELL", "CD27",
                    "GZMB", "PRF1", "TBX21", "PRDM1", "NKG7", "GNLY")
myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
library(ggrepel)
pdf(file = "../Figures/Figure_2/Figure2B.pdf", width = 6, height = 6)
ggplot(data = myData, aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(alpha = 0.2, size = 1) +
  theme(legend.position = "none", text = element_text(size = 10)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ## Colours for significance
  geom_point(data = subset(myData, padjThresh==T & logFC<(-0.5)), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "blue4") +
  geom_point(data = subset(myData, padjThresh==T & logFC>0.5), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "orangered3") +
  ## Label these
  geom_label_repel(mapping = aes(label = Labels), box.padding = 0.35, point.padding = 0.2) +
  # geom_text_repel(data=subset(myData,padjThresh==TRUE & logFC>0.5), aes(logFC,-log10(P.Value),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
  labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
dev.off()

## Quick check to see that numbers match those in Figure 2C
# nrow(subset(myData, padjThresh==T & logFC<(-0.5))) + nrow(subset(myData, padjThresh==T & logFC>0.5))




# Figure 2C - Venn diagram of differentially expressed genes between VD1 and CD8 compartments
de.common <- which(dt[, "CD27LOvsCD27HI"]!=0 & dt[,"EMRAvsNaive"]!=0)
de.VD1 <- which(dt[, "CD27LOvsCD27HI"]!=0 & dt[, "EMRAvsNaive"]==0)
de.CD8 <- which(dt[, "CD27LOvsCD27HI"]==0 & dt[, "EMRAvsNaive"]!=0)

length(de.common) # Input these values into http link below (uncomment and Cmd + click)
length(de.VD1)
length(de.CD8)
# http://eulerr.co/
Common_genes <- tfit$genes$ENTREZ[de.common]


# Figure 2D - Heatmap of those genes which are DE in common between populations
## Remove VD2
v1 <- v
check <- as.data.frame(v$E)
v1$targets <- droplevels(subset(v1$targets, group != "VD2"))
v1$E <- as.data.frame(v1$E)
v1$E <- v1$E %>% dplyr::select(-contains("VD2"))
v1$E <- as.matrix(v1$E)
j <- !grepl("VD2", group)
group1 <- group[j]

## subset to only contain genes in Common_genes
i <- which(v1$genes$ENTREZID %in% Common_genes)

## Define needed bits
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999")[v1$targets$group]
data.frame(v1$targets$group, col.cell1)
mycol <- colorpanel(1000,"blue","white","red")

pdf("../Figures/Figure_2/Figure2D_vertical_P_Yl.pdf", width = 10, height = 50)
heatmap.2(v1$E[i,], scale = "row",
          labRow = ifelse((is.na(v1$genes$SYMBOL[i])), v1$genes$ENTREZID[i], v1$genes$SYMBOL[i]), labCol = NA,
          col = Seurat::PurpleAndYellow(), trace = "none", density.info = "none", 
          margin = c(8, 6), lhei = c(2, 10), 
          hclustfun = hclustAvg, ColSideColors = col.cell1)
dev.off()


# Figure 2E - Cytoscape of pathways
## Open connections
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


library(qusage)

## Read in the Genesets
All_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/msigdb.v6.2.entrez.gmt")
KEGG_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/c2.cp.kegg.v6.2.entrez.gmt")
GO_terms <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Output/Cyto/c5.all.v6.2.entrez.gmt")


## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_terms)]
filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_terms)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (4326 genesets)
geneset_sizes <- unlist(lapply(filter_gmt, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- filter_gmt[geneset_indices]

## Perform camera analysis on filtered geneset
idx <- ids2indices(filtered_set, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
nrow(camera_results)
droplevels(subset(camera_results, FDR <= 0.001))


# BiocManager::install("qusage")
library(qusage)
camera_results_a <- camera_results[rownames(camera_results) %in% names(filtered_set),]

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

camera_results_generic_em <- data.frame(rownames(camera_results_a), camera_descr, 
                                        PValue = camera_results_a[, "PValue"],
                                        FDR = camera_results_a[, "FDR"],
                                        camera_Phenotype,
                                        camera_genes)

camera_results_file <- "../../Cyto/camera_results_generic_CD8.txt"
write.table(camera_results_generic_em, file.path(camera_results_file), 
            col.name = T, sep = "\t", row.names = F, quote = F)
# 
# camera_results_CD8 <- droplevels(subset(camera_results_generic_em, FDR <= 0.001))
# camera_results_CD8$camera_Phenotype <- ifelse((camera_results_CD8$camera_Phenotype == 1), "CD8 EMRA", "CD8 Naive")
# write.csv("../../Cyto/CD8_results.csv", x = camera_results_CD8, row.names = F)
# 
# camera_results_VD1$camera_Phenotype <- ifelse((camera_results_VD1$camera_Phenotype == 1), "VD1 CD27LO", "VD1 CD27HI")
# write.csv("../../Cyto/VD1_results.csv", x = camera_results_VD1, row.names = F)

expression_file <- "../../Cyto/expression_file.txt"
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
analysis_name <- "EMRA_vs_Naive"
cur_model_name <- paste("camera", analysis_name, sep="_")
results_filename <- file.path(getwd(),  camera_results_file)


current_network_name <- paste(cur_model_name, pvalue_threshold, qvalue_threshold, sep = "_")


results_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Cyto/camera_results_generic_CD8.txt"
expression_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Cyto/expression_file.txt"
generic_gmt_file <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Cyto/msigdb.v6.2.entrez.gmt"
em_command = paste('enrichmentmap build analysisType=generic',
                   "gmtFile=", generic_gmt_file,
                   "pvalue=", pvalue_threshold,
                   "qvalue=", qvalue_threshold,
                   "similaritycutoff=", similarity_threshold,
                   "coefficients=", similarity_metric,
                   "enrichmentsDataset1=", results_filename,
                   "expressionDataset1=", expression_filename)

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


###### VD2 Stuff
## Hierarchical clustering - no longer produces nice one
DF <- as.data.frame(lcpm)

DF1 <- #DF[, !grepl( "VD2" , names(DF) )] 
  DF %>% 
  rownames_to_column(var = "Genes") %>% 
  gather(contains("."), key = "ID" , value = "value") %>%
  spread(key = "Genes", value = "value") 

DF2 <- data.frame(DF1[, names(DF1) != "ID"], row.names = DF1[, names(DF1) == "ID"])
hc <- hclust(dist(DF2), "average")
p1 <- ggdendrogram(hc, rotate = F, size = 2, leaf_labels = F)

df2 <- data.frame(cluster = cutree(hc, 2), cell.types = factor(hc$labels, levels = hc$labels[hc$order]))

merging <- as.data.frame(x$samples) 
merging1 <- rownames_to_column(merging, var = "cell.types")

merging2 <- merging1[, c("cell.types", "group")]

df3 <- merge(df2, merging2, by = "cell.types")

p2 <- ggplot(df3, aes(cell.types, y = 1, fill = factor(cluster))) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), 
                    name = "Cluster",
                    breaks = c("1", "2"),
                    labels = c("Cluster 1", "Cluster 2"))
p3 <- ggplot(df3, aes(cell.types, y = 1, fill = factor(group))) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = cbcols, 
                    name = "Cluster")
p3
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)  
gp3 <- ggplotGrob(p3)

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

library(gridExtra)
pdf("../../Bulk/Figures/VD2/Dendro.pdf")
grid.arrange(gp1, gp3, ncol = 1, heights = c(4/5, 1/5, 1/5)) # Doesn't 
dev.off()


col.cell <- c("#56B4E9","#E69F00","#009E73","#999999", "#CC79A7")[x$samples$group]
data.frame(x$samples$group, col.cell)

pdf("../Figures/VD2/MDS_scaling.pdf", height = 6, width = 6)
plotMDS(lcpm, pch = 16, cex = 2, col = col.cell, top = 14191)
legend("top",
       fill = c("#56B4E9",
                "#E69F00", "#009E73", "#999999", "#CC79A7"),
       legend = levels(x$samples$group))
dev.off()



## subset to only contain genes in Common_genes
de.common <- which(dt[, "CD27LOvsCD27HI"]!=0 & dt[,"EMRAvsNaive"]!=0)
Common_genes <- tfit$genes$ENTREZID[de.common]
i <- which(v$genes$ENTREZID %in% Common_genes)

## Define needed bits
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")[v$targets$group]
data.frame(v$targets$group, col.cell1)
mycol <- colorpanel(1000,"blue","white","red")

pdf("../Figures/VD2/426_genes_VD2.pdf", width = 10, height = 50)
heatmap.2(v$E[i,], scale = "row",
          labRow = ifelse((is.na(v$genes$SYMBOL[i])), v$genes$ENTREZID[i], v$genes$SYMBOL[i]), labCol = NA,
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8, 6), lhei = c(2, 10), 
          hclustfun = hclustAvg, ColSideColors = col.cell1)
dev.off()


VD2_spe <- which(dt[, "CD27HIvsVD2"]!=0 & dt[, "EMRAvsVD2"]!=0 & 
                   dt[, "CD27LOvsVD2"]!=0 & dt[, "NaivevsVD2"]!=0)
tfit$genes$SYMBOL[VD2_spe]

effector.common <- which(dt[, "CD27HIvsVD2"]!=0 & dt[,"CD27LOvsVD2"]==0)
Common_genes <- tfit$genes$SYMBOL[effector.common]
Common_genes
naive.common <- which(dt[, "CD27HIvsVD2"]==0 & dt[,"CD27LOvsVD2"]!=0)
n_genes <- tfit$genes$SYMBOL[naive.common]
n_genes


VD2_eff <- droplevels(CD27HI.vs.VD2[CD27HI.vs.VD2$SYMBOL %in% Common_genes,])


CD27HI.vs.VD2[grepl("GZMA", CD27HI.vs.VD2$SYMBOL),]
CD27LO.vs.VD2[grepl("CCR2", CD27LO.vs.VD2$SYMBOL),]


######## OTher Stuff #####

### GSEA DE genes

barcodeplot(tfit$t[, "CD27LOvsCD27HI"], index = idx$GO_CELL_KILLING, main = ".")
barcodeplot(tfit$t[, "CD27LOvsCD27HI"], index = idx$KEGG_OXIDATIVE_PHOSPHORYLATION, main = ".")
barcodeplot(tfit$t[, "CD27LOvsCD27HI"], index = idx, main = ".")








## PCA on the 'effector population' for all genes... What does PC1 correspond to?
### Remove Naive populations
v2 <- v
v2$targets <- droplevels(subset(v2$targets, group != "VD1.CD27HI" & group != "CD8.Naive"))
v2$E <- as.data.frame(v2$E)
v2$E <- v2$E %>% dplyr::select(-contains("VD1.CD27HI")) %>% dplyr::select(-contains("CD8.Naive"))
v2$E <- as.matrix(v2$E)

## Do a grepl on multiple things...
patterns <- c("CD8.Naive", "VD1.CD27HI")
j <- !grepl(paste(patterns, collapse="|"), group)
group2 <- group[j]
i <- v2$genes$ENTREZID


this <- as.data.frame(v2$E[i,])
that <- rownames_to_column(this, var = "Gene")
that1 <- that %>% gather(contains("L00"), key = "Sample", value = "value")
that2 <- spread(that1, key = Gene, value = value)

# Label groups
subgroup <- ifelse(grepl("VD1.CD27LO", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", that2$Sample), "VD1.CD27LO",
                   ifelse(grepl("CD8.EMRA", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", that2$Sample), "CD8.EMRA",
                          ifelse(grepl("CD8.Naive", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", that2$Sample), "CD8.Naive",
                                 ifelse(grepl("VD1.CD27HI", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", that2$Sample), "VD1.CD27HI",
                                        ifelse(grepl("VD2", that2$Sample), "VD2", "none")))))
that2$subgroup <- as.factor(subgroup)

# Remove identifiers
that3 <- data.frame(that2[, names(that2) != "Sample"], row.names = that2[, names(that2) == "Sample"])
colnames(that3) <- gsub("^X", "",  colnames(that3))

prin_comp <- prcomp(that3[, names(that3) != "subgroup"], scale. = T)
subgroup <- that3[, "subgroup"]


g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = subgroup, ellipse = F,
              circle = T,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = "horizontal", 
               legend.position = "top")

g <- g + ggtitle("PCA of shared DE genes (Naive vs Effector)")




##### SUPPLEMENTARY FIGURES ####
Common_genes <- tfit$genes$ENTREZ[de.common]
VD1_genes <- tfit$genes$ENTREZ[de.VD1]
CD8_genes <- tfit$genes$ENTREZ[de.CD8]

int <- as.vector(matrix(c(Common_genes, VD1_genes), nrow = 2, byrow = T)) 
all_gen <- as.vector(matrix(c(int, CD8_genes), nrow = 2, byrow = T)) 

allgen <- all_gen[!duplicated(all_gen)]

# Venn Diagrams - S2a
png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/v2/DE_genes_EffectorvsNaive.png",
    width = 300, height = 300, units = "mm", res = 300)
vennDiagram(dt[, c(1,8)], include = "both", mar = c(0,0,0,0))
title(main = "Differentially expressed Genes", line = -6)
dev.off()

png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/v2/Up_DE_genes_EffectorvsNaive.png",
    width = 300, height = 300, units = "mm", res = 300)
vennDiagram(dt[, c(1,8)], include = "up", mar = c(0,0,0,0))
title(main = "Upregulated", line = -6)
dev.off()

png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/v2/Down_DE_genes_EffectorvsNaive.png",
    width = 300, height = 300, units = "mm", res = 300)
vennDiagram(dt[, c(1,8)],  include = "down", mar = c(0,0,0,0))
title(main = "Downregulated", line = -6)
dev.off()








## Get data in correct format
de.Eff <- which(dt[, "CD27LOvsEMRA"]!=0)
Symbol <- tfit$genes$SYMBOL[de.Eff]

de.Naive <- which(dt[, "CD27HIvsNaive"]!=0)
Symbol <- tfit$genes$SYMBOL[de.Naive]
length(Symbol)

View(Symbol)

Symbol1 <- Symbol[!grepl("^TRAV", Symbol)]
Symbol1 <- Symbol1[!grepl("^TRBV", Symbol1)]
Symbol1 <- Symbol1[!grepl("^TRD", Symbol1)]
Symbol1 <- Symbol1[!grepl("^TRG", Symbol1)]
Symbol1 <- Symbol1[!grepl("^TRAC", Symbol1)]


camera_results_naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27HIvsNaive"])
head(camera_results_naive)

sig.cam <- subset(camera_results_naive, PValue < 0.05 & FDR < 0.0005)
sig.cam

myData <- as.data.frame(CD27HI.vs.Naive)
myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)
myData$Labels <- myData$SYMBOL
myData$Labels <- ""
labelled_genes <- c("LTB", "LEF1", "TCF7", "CCR7", "SELL", "CD27",
                    "GZMB", "PRF1", "TBX21", "PRDM1", "NKG7", "GNLY")
myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""

# pdf(file = "../Figures/Figure_2/Figure2B.pdf", width = 6, height = 6)
ggplot(data = myData, aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(alpha = 0.2, size = 1) +
  theme(legend.position = "none", text = element_text(size = 10)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ## Colours for significance
  geom_point(data = subset(myData, padjThresh==T & logFC<(-0.5)), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "blue4") +
  geom_point(data = subset(myData, padjThresh==T & logFC>0.5), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "orangered3") +
  ## Label these
  geom_label_repel(mapping = aes(label = Labels), box.padding = 0.35, point.padding = 0.2) +
  # geom_text_repel(data=subset(myData,padjThresh==TRUE & logFC>0.5), aes(logFC,-log10(P.Value),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
  labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
dev.off()


CD27LO.vs.VD2$Direction <- ifelse((CD27LO.vs.VD2$logFC <= 0), "VD2", "CD27LO")
sig <- droplevels(subset(CD27LO.vs.VD2, adj.P.Val <= 0.05))
CD27LO.vs.VD2[grepl("TRDV2", CD27LO.vs.VD2$SYMBOL), ]
write.csv("../Output/DEgenes/CD27LO_vs_VD2.csv", x = CD27LO.vs.VD2, row.names = F)

CD27HI.vs.VD2$Direction <- ifelse((CD27HI.vs.VD2$logFC <= 0), "VD2", "CD27HI")

write.csv("../Output/DEgenes/CD27HI_vs_VD2.csv", x = CD27HI.vs.VD2, row.names = F)
sig <- droplevels(subset(CD27HI.vs.VD2, adj.P.Val <= 0.05))





CD27HI.vs.Naive$Direction <- ifelse((CD27HI.vs.Naive$logFC <= 0), "CD8_Naive", "CD27HI")



head(CD27HI.vs.Naive)
write.csv("../Output/DEgenes/CD27HI.vs.Naive.csv", x = CD27HI.vs.Naive, row.names = F)




head(EMRA.vs.Naive)

this <- droplevels(subset(EMRA.vs.Naive, adj.P.Val < 0.05))
head(this)

nrow(this)
length(Common_genes)

try <- this[!'%in%'(this$ENTREZID, Common_genes), ]
head(try)
nrow(try)

try$Direction <- ifelse((try$logFC <= 0), "CD8_Naive", "CD8_EMRA")

nrow(droplevels(subset(try, Direction == "CD8_Naive")))

658/1328



