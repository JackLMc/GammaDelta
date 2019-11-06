library(tidyverse)
## Looking at how many genes are Transcription factors
VD1 <- read.csv("./Bulk/Output/CD27LO.vs.CD27HI.csv")
VD1a <- droplevels(subset(VD1, adj.P.Val <= 0.05))

TF <- read.csv("./Transcription_Factor/TFCheckpoint.csv")
colnames(TF)[colnames(TF) == "gene_symbol"] <- "SYMBOL"
colnames(TF)[colnames(TF) == "entrez_human"] <- "ENTREZID"
TF_imp <- droplevels(subset(TF, AnimalTFDB == "AnimalTFDB"))
TF_imp <- TF_imp[, c("SYMBOL", "ENTREZID", "gene_name", "AnimalTFDB")]
VD1b <- merge(VD1a, TF_imp[, c("SYMBOL", "ENTREZID", "gene_name")],
                by = c("SYMBOL", "ENTREZID")) %>% droplevels()

VD1b$Direction <- ifelse((VD1b$logFC > 0), "VD1.CD27LO", "VD1.CD27HI")
write.csv("./Transcription_Factor/VD1_DGE_TF.csv", x = VD1b, row.names = F)

### CD8
CD8 <- read.csv("./Bulk/Output/EMRA.vs.Naive.csv")
CD8a <- droplevels(subset(CD8, adj.P.Val <= 0.05))

TF <- read.csv("./Transcription_Factor/TFCheckpoint.csv")
colnames(TF)[colnames(TF) == "gene_symbol"] <- "SYMBOL"
colnames(TF)[colnames(TF) == "entrez_human"] <- "ENTREZID"
TF_imp <- droplevels(subset(TF, AnimalTFDB == "AnimalTFDB"))
TF_imp <- TF_imp[, c("SYMBOL", "ENTREZID", "gene_name", "AnimalTFDB")]
CD8b <- merge(CD8a, TF_imp[, c("SYMBOL", "ENTREZID", "gene_name")],
              by = c("SYMBOL", "ENTREZID")) %>% droplevels()

CD8b$Direction <- ifelse((CD8b$logFC > 0), "CD8.EMRA", "CD8.Naive")
write.csv("./Transcription_Factor/CD8_DGE_TF.csv", x = CD8b, row.names = F)

load("./Bulk/Final_edgeR.RData")

## Versus VD2
TF_imp <- droplevels(subset(TF, AnimalTFDB == "AnimalTFDB"))
head(TF_imp)
sig <- droplevels(subset(CD27LO.vs.VD2, adj.P.Val <= 0.05))
sig$SYMBOL <- as.factor(sig$SYMBOL)
TF_imp$SYMBOL <- as.factor(TF_imp$SYMBOL)
tf_VD2_CD27LO <- droplevels(sig$SYMBOL[sig$SYMBOL %in% TF_imp$SYMBOL])

sig <- droplevels(subset(CD27HI.vs.VD2, adj.P.Val <= 0.05))
sig$SYMBOL <- as.factor(sig$SYMBOL)
TF_imp$SYMBOL <- as.factor(TF_imp$SYMBOL)
tf_VD2_CD27HI <- droplevels(sig$SYMBOL[sig$SYMBOL %in% TF_imp$SYMBOL])

write.csv("./Transcription_Factor/VD2_CD27LO.csv", x = tf_VD2_CD27LO)
write.csv("./Transcription_Factor/VD2_CD27HI.csv", x = tf_VD2_CD27HI)


# Comparisons and Colours
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"), c("VD1.CD27HI", "CD8.EMRA"), c("VD1.CD27HI", "CD8.Naive"), c("VD1.CD27HI", "VD2"),
                       c("VD1.CD27LO", "CD8.EMRA"), c("VD1.CD27LO", "CD8.Naive"), c("VD1.CD27LO", "VD2"),
                       c("CD8.EMRA", "CD8.Naive"), c("CD8.EMRA", "VD2"), c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999", "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00", "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")

## The transcription factor heatmap and bespoke
## Do a heatmap based on the VD1 transcription factors
library(gplots)
mycol <- colorpanel(1000,"blue","white","red")
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

## v has VD2, v1 has not
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")[v1$targets$group]
i <- which(v$genes$SYMBOL %in% VD1b$SYMBOL)
pdf(file = "./Transcription_Factor/VD1_TF_heatmap.pdf", height = 8, width = 8)
heatmap.2(v1$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = NA,
          col = mycol,
          trace = "none", density.info = "none",
          margin = c(10,6), lhei = c(2,10), ColSideColors = col.cell1)
dev.off()



