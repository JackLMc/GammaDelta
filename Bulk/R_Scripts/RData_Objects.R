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


load("Bulk/Final_edgeR.RData")

# Targeted
use_these <- c("CD28", "CD27", "JAML", "TNFSF8", 
               "KLRK1", "NCR1", "KLRF1",
               "TIGIT", "SLAMF7", "PDCD1", "LAG3")

# NKG2D: KLRK1, NKP46: NCR1, NKP80: KLRF1
# mat_data <- t(v1$E[i,])
# conversions <- droplevels(v1$genes[v1$genes$ENTREZID %in% colnames(mat_data), ])

i <- which(v1$genes$SYMBOL %in% use_these)
pdf(file = "./Paper/targeted_heatmap.pdf", height = 8, width = 8)
heatmap.2(t(v1$E[i,]), scale = "column",
          labCol = v$genes$SYMBOL[i], labRow = NA,
          col = mycol, dendrogram = "row", 
          trace = "none", density.info = "none",
          margin = c(10,6), lhei = c(2,10),
          hclustfun = hclustAvg, 
          RowSideColors = col.cell1)
dev.off()

# Making RData objects of genesets# Soruce script
# source("./Bulk/R_Scripts/Counts.R")


## JUST NEED A LIST NOW!
Carrie <- read.csv("/Volumes/ResearchData/Willcox Group/Jack/gamma_delta_RNASeq/Genesets_to_investigate/Carrie_Geneset 150519 with Vd2.csv")
head(Carrie)
Carrie$GENESET <- as.factor(Carrie$GENESET)
Carrie_RData <- list()
c <- 1
for(i in levels(Carrie$GENESET)){
  print(i)
  work <- droplevels(subset(Carrie, GENESET == i))
  vec <- as.vector(work$SYMBOL)
  Carrie_RData[[i]] <- vec
  c <- c + 1
}

# Heatmaps
## source("Counts.R")
library(gplots)
mycol <- colorpanel(1000,"blue","white","red")
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")[v$targets$group]
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999")[v1$targets$group]

# c <- 1
# for(j in levels(Carrie$GENESET)){
#   print(j)
#   working <- droplevels(subset(Carrie, GENESET == j))
#   i <- which(v1$genes$SYMBOL %in% working$SYMBOL)
#   setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Heatmaps/V4")
#   filen <- paste0(j, ".pdf")
#   pdf(filename = filen, height = 8, width = 8)
# heatmap.2(v1$E[i,], scale = "row",
#           labRow = v1$genes$SYMBOL[i], labCol = colnames(v1$E),
#           col = mycol,
#           trace = "none", density.info = "none",
#           margin = c(10,6), lhei = c(2,10),
#           #dendrogram = "column",
#           main = as.character(levels(working$GENESET)),
#           hclustfun = hclustAvg)
# dev.off()
#   setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk")
#   c <- c + 1
# }


Carrie <- read.csv("/Volumes/ResearchData/Willcox Group/Jack/gamma_delta_RNASeq/Genesets_to_investigate/Carrie_Geneset 150519 with Vd2.csv")
head(Carrie)
Carrie$GENESET <- as.factor(Carrie$GENESET)
Carrie_RData <- list()
c <- 1
for(i in levels(Carrie$GENESET)){
  print(i)
  work <- droplevels(subset(Carrie, GENESET == i))
  vec <- as.vector(work$SYMBOL)
  Carrie_RData[[i]] <- vec
  c <- c + 1
}

c <- 1
for(j in levels(Carrie$GENESET)){
  print(j)
  working <- droplevels(subset(Carrie, GENESET == j))
  i <- which(v$genes$SYMBOL %in% working$SYMBOL)
  setwd("/Volumes/ResearchData/Willcox Group/Jack/gamma_delta_RNASeq/Figures/Heatmaps/Merged_Lanes/V5/With_VD2/")
  filen <- paste0(j, ".pdf")
  pdf(file = filen, height = 8, width = 8)
  heatmap.2(v$E[i,], scale = "row",
            labRow = v$genes$SYMBOL[i], labCol = NA,
            col = mycol,
            trace = "none", density.info = "none",
            margin = c(10,6), lhei = c(2,10),
            #dendrogram = "column",
            main = as.character(levels(working$GENESET)),
            ColSideColors = col.cell)
  dev.off()
  setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk")
  c <- c + 1
}



 save.image("Carrie_Heatmap.RData")


# Removed VD2
Carrie <- read.csv("/Volumes/ResearchData/Willcox Group/Jack/gamma_delta_RNASeq/Genesets_to_investigate/Carrie_Geneset 150519 without Vd2.csv")
head(Carrie)
Carrie$GENESET <- as.factor(Carrie$GENESET)
Carrie_RData <- list()
c <- 1
for(i in levels(Carrie$GENESET)){
  print(i)
  work <- droplevels(subset(Carrie, GENESET == i))
  vec <- as.vector(work$SYMBOL)
  Carrie_RData[[i]] <- vec
  c <- c + 1
}


c <- 1
for(j in levels(Carrie$GENESET)){
  print(j)
  working <- droplevels(subset(Carrie, GENESET == j))
  i <- which(v1$genes$SYMBOL %in% working$SYMBOL)
  setwd("/Volumes/ResearchData/Willcox Group/Jack/gamma_delta_RNASeq/Figures/Heatmaps/Merged_Lanes/V4/Without_VD2")
  filen <- paste0(j, ".pdf")
  pdf(file = filen, height = 8, width = 8)
  heatmap.2(v1$E[i,], scale = "row",
            labRow = v1$genes$SYMBOL[i], labCol = NA, 
            col = mycol,
            trace = "none", density.info = "none", 
            margin = c(10,6), lhei = c(2,10), 
            #dendrogram = "column",
            main = as.character(levels(working$GENESET)),
            ColSideColors = col.cell1)
  dev.off()
  setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk")
  c <- c + 1
}






dev.off()

i <- which(v1$genes$SYMBOL %in% Carrie$SYMBOL)

mycol <- colorpanel(1000,"blue","white","red")

included <- rownames(v1$E[i, ])

these <- v1$genes[v1$genes$ENTREZID %in% included, ]
Carrie <- Carrie[!duplicated(Carrie$SYMBOL), ]

that <- merge(these, Carrie, by = "SYMBOL")
col.cell1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[that$GENESET]
data.frame(that$GENESET, col.cell1)

heatmap.2(v$E[i,], scale = "row",
          labRow = v1$genes$SYMBOL[i], labCol = colnames(v$E), 
          col = mycol, RowSideColors = col.cell1,
          trace = "none", density.info = "none", 
          margin = c(10,6), lhei = c(2,10), 
          #dendrogram = "column",
          hclustfun = hclustAvg)



