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




# Bulk RNA Sequencing Analysis of Gamma delta cell #

## Introduction ##
#### This is the final set of transcriptomic analyses of the gamma-delta cells. 
#### This utilises the sorted BAM files which were produced using HISAT2
#### This also utilises the counts created by the CreateCounts.R script
#### CreateCounts.R script utilses [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)

## Loading tools and data ##
#### We're using the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package for analysis
source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR", dependencies = T)
# biocLite("Homo.sapiens", dependencies = T)
source("Bulk/R_Scripts/Functions.R")
# source("CreateCounts.R")

## Using edgeR ##
library(edgeR)
files <- list.files(path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Counts", pattern = ".txt$")
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Counts")

#### Create the DGE List object used to investigate
x <- readDGE(files, columns = c(1, 3))
samplenames <- substring(colnames(x), length(files), nchar(colnames(x)))
colnames(x) <- samplenames

## Preprocessing of data ##
#### Label the data based on cell type and lane.
#### Rename the columns from the extra long names provided, to shorter names
group <- ifelse(grepl("VD1.CD27LO", colnames(x)) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", colnames(x)), "VD1.CD27LO",
                ifelse(grepl("CD8.EMRA", colnames(x)) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", colnames(x)), "CD8.EMRA",
                       ifelse(grepl("CD8.Naive", colnames(x)) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", colnames(x)), "CD8.Naive",
                              ifelse(grepl("VD1.CD27HI", colnames(x)) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", colnames(x)), "VD1.CD27HI",
                                     ifelse(grepl("VD2", colnames(x)), "VD2", "none")))))
x$samples$group <- as.factor(group)

lane <- ifelse(grepl("L001", colnames(x)), "L001",
               ifelse(grepl("L002", colnames(x)), "L002",
                      ifelse(grepl("L003", colnames(x)), "L003",
                             ifelse(grepl("L004", colnames(x)), "L004", "none"))))
x$samples$lane <- as.factor(lane)

colnames(x) <- c("28MD.VD1.CD27LO.L001", "28MD.VD1.CD27LO.L002", "28MD.VD1.CD27LO.L003", "28MD.VD1.CD27LO.L004",
                 "28MD.CD8.EMRA.L001", "28MD.CD8.EMRA.L002","28MD.CD8.EMRA.L003", "28MD.CD8.EMRA.L004",
                 "28MD.CD8.Naive.L001", "28MD.CD8.Naive.L002", "28MD.CD8.Naive.L003", "28MD.CD8.Naive.L004",
                 "28MD.VD1.CD27HI.L001", "28MD.VD1.CD27HI.L002", "28MD.VD1.CD27HI.L003", "28MD.VD1.CD27HI.L004",
                 "28MD.VD2.L001", "28MD.VD2.L002", "28MD.VD2.L003", "28MD.VD2.L004",
                 "31EN.CD8.Naive.L001", "31EN.CD8.Naive.L002", "31EN.CD8.Naive.L003", "31EN.CD8.Naive.L004",
                 "31EN.CD8.EMRA.L001", "31EN.CD8.EMRA.L002", "31EN.CD8.EMRA.L003", "31EN.CD8.EMRA.L004",
                 "31EN.VD1.CD27LO.L001", "31EN.VD1.CD27LO.L002", "31EN.VD1.CD27LO.L003", "31EN.VD1.CD27LO.L004",
                 "31EN.VD1.CD27HI.L001", "31EN.VD1.CD27HI.L002", "31EN.VD1.CD27HI.L003", "31EN.VD1.CD27HI.L004",
                 "31EN.VD2.L001", "31EN.VD2.L002", "31EN.VD2.L003", "31EN.VD2.L004")
samplenames <- colnames(x)

#### Rename the gene names from entrez identifiers to gene symbols to make it easier to interpret results
#### Some genes map to more than one Entrez ID, remove those which are duplicated
library(Homo.sapiens)
geneid <- rownames(x)
genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM"), 
                keytype = "ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes


#### We now have a DGEList with gene names attached which are symbols
#### Calculate counts per million and log counts per million
#### This normalises for library sizes and filters out lowly expressed transcripts
cpm_ <- cpm(x$counts)

lcpm <- cpm(x$counts, log = T)
View(cpm_)

#### Some genes are very lowly expressed in a lot of samples, remove these to reduce computing power later
#### Show me the amount of transcripts that are zero for all 40 samples
table(rowSums(x$counts==0)==40)

#### Genes must have a cpm above 1.8 (count of 6.5 in lowest library)
#### and be expressed in at least 8 runs to be kept (one group)
keep.exprs <- rowSums(cpm_>1.8)>=8 
x <- x[keep.exprs,, keep.lib.sizes = F]

#### There are ~12,000 genes which pass these criteria
dim(x)

#### Normalisation for RNA composition is required
#### Avoid over-estimation of differentially expressed genes
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#### There are no batch effects in the data, and therefore this has not been completed
#### all.equal(removeBatchEffect(x), x$counts)

#### Estimating dispersion
d1 <- estimateCommonDisp(x, verbose = T)
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + x$samples$group)
colnames(design.mat) <- levels(x$samples$group)
d2 <- estimateGLMCommonDisp(x,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method = "power")

d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

save.image("/Final_edgeR.RData")
load("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Final_edgeR.RData")

#### Investigate similarity of cell populations
library(ggdendro) 
library(tidyverse)
DF <- as.data.frame(lcpm)

DF1 <- #DF[, !grepl( "VD2" , names(DF) )] 
  DF %>% 
  rownames_to_column(var = "Genes") %>% 
  gather(contains("."), key = "ID" , value = "value") %>%
  spread(key = "Genes", value = "value") 

DF2 <- data.frame(DF1[, names(DF1) != "ID"], row.names = DF1[, names(DF1) == "ID"])
hc <- hclust(dist(DF2), "ave")
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

pdf("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/hier_clust_VD2.pdf")
grid.arrange(gp1, gp3, ncol = 1, heights = c(4/5, 1/5, 1/5))
dev.off()

#### Multidimensional scaling
No_VD2 <- droplevels(x$samples$group[x$samples$group != "VD2"])
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999")[No_VD2]
data.frame(No_VD2, col.cell1)
DFa <- DF[, !grepl( "VD2" , names(DF) )]

pdf("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Figure2A.pdf")
plotMDS(DF, pch = 16, cex = 1.2, col = cbcols)
legend("topleft",
       fill = c("#56B4E9", "#E69F00",
                "#009E73", "#999999"),
       legend = levels(No_VD2))
dev.off()

