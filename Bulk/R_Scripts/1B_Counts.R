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


library(BiocManager)
# BiocManager::install("Rsubread", dependencies = T)
library(Rsubread)
## If not previously ran, CreateCounts by sourcing...
# source("CreateCounts.R") 

# use edgeR  START!
# BiocManager::install("edgeR", dependencies = T)
library(edgeR)
y <- read.delim("./Bulk/Counts/Counts_GD.txt", check.names = F) %>% column_to_rownames(., var = "GeneID") %>% as.matrix

# Replenish "x"
x <- DGEList(counts = y)

# Need to label based on cell type and lane
## Cell types
colnames(x) # BP28-VD1 is CD27LO
# group <- ifelse(grepl("VD1.CD27LO", colnames(x), ignore.case = T) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", colnames(x), ignore.case = T), "VD1.CD27LO",
#                           ifelse(grepl("CD8.EMRA", colnames(x), ignore.case = T) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", colnames(x), ignore.case = T), "CD8.EMRA",
#                                  ifelse(grepl("CD8.Naive", colnames(x), ignore.case = T) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", colnames(x), ignore.case = T), "CD8.Naive",
#                                         ifelse(grepl("VD1.CD27HI", colnames(x), ignore.case = T) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", colnames(x), ignore.case = T), "VD1.CD27HI",
#                                                ifelse(grepl("VD2", colnames(x), ignore.case = T), "VD2", "none")))))

colnames(x) <- gsub("3_BP28Vd1.bam", "3_BP28Vd1CD27LO.bam", colnames(x))

group <- ifelse(grepl("CD27LO", colnames(x), ignore.case = T), "VD1.CD27LO",
                          ifelse(grepl("EMRA", colnames(x), ignore.case = T) , "CD8.EMRA",
                                 ifelse(grepl("Naive", colnames(x), ignore.case = T), "CD8.Naive",
                                        ifelse(grepl("CD27HI", colnames(x), ignore.case = T), "VD1.CD27HI",
                                               ifelse(grepl("VD2", colnames(x), ignore.case = T), "VD2", "none")))))

x$samples$group <- as.factor(group)
colnames(x) <- c("BP19.CD8.EMRA", "BP26.VD1.CD27LO", "BP29.VD1.CD27LO", "CD22.CD8.EMRA", "BP29.CD8.Naive", "BP29.CD8.EMRA",
                 "BP27.CD8.Naive", "BP27.CD8.EMRA", "BP27.VD1.CD27HI", "BP17.VD1.CD27LO", "BP19.CD8.Naive", "BP26.CD8.EMRA",
                 "28MD.CD8.EMRA", "28MD.CD8.Naive", "28MD.VD1.CD27HI", "28MD.VD1.CD27LO", "28MD.VD2", "BP24.CD8.Naive", "BP28.VD1.CD27LO",
                 "BP9.VD1.CD27HI", "31EN.CD8.EMRA", "31EN.CD8.Naive", "31EN.VD1.CD27HI", "31EN.VD1.CD27LO", "31EN.VD2", "BP24.CD8.EMRA",
                 "BP29.VD1.CD27HI", "BP28.CD8.Naive", "BP28.CD8.EMRA", "BP22.CD8.Naive")
samplenames <- colnames(x)

# Gaining second dataframe (Symbols)
# BiocManager::install("Homo.sapiens",dependencies = T)
library(Homo.sapiens)
geneid <- rownames(x)
genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM"), 
                keytype = "ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes


# Calculate counts per million, log counts per million (can also do RPKM [function rpkm()])
cpm <- cpm(x)
lcpm <- cpm(x, log = T)

# Remove lowly expressed transcripts
dim(x)

unique(samplenames) %>% length()
table(rowSums(x$counts==0)==30) # 5,377 genes are zero in all 30 samples


CPM_scaling <- min(x$samples$lib.size)/1000000
cpm_scale <- 6.5/CPM_scaling

y <- x



keep.exprs <- rowSums(cpm>cpm_scale)>=2 # Genes must have a cpm above 0.44 (count of 6.5 in lowest library) and be expressed in at least 2 groups (1 population)
x <- x[keep.exprs,, keep.lib.sizes = F]
dim(x) # 12,616 kept


# What genes have I lost
pre <- y@.Data[3] %>% as.data.frame()
post <- x@.Data[3] %>% as.data.frame()

filtered_genes <- pre[!('%in%'(pre$ENTREZID, post$ENTREZID)), ] %>% droplevels() # 15,779 filtered


# ## Showing the removal of the lowly expressed transcripts
# library(RColorBrewer)
# nsamples <- ncol(x)
# col <- brewer.pal(nsamples)
# pdf("./Bulk/Figures/Proof/lowly_expressed_transcripts.pdf")
# par(mfrow = c(1,2))
# plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2,
#      main = "", xlab = "")
# title(main = "A. Raw data", xlab = "Log-cpm")
# abline(v = 0, lty = 3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col = col[i], lwd = 2)
# }
# legend("topright", samplenames, text.col = col, bty = "n")
# lcpm <- cpm(x, log = T)
# plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2,
#      main = "", xlab = "")
# title(main = "B. Filtered data", xlab = "Log-cpm")
# abline(v = 0, lty = 3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col = col[i], lwd = 2)
# }
# legend("topright", samplenames, text.col = col, bty = "n")
# dev.off()

# Normalising gene expression
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
x <- estimateCommonDisp(x)
x <- estimateTagwiseDisp(x)
# save(x, file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/NormalisationData.RData")

# ## Showing theoretical effect of normalisation
# x2 <- x
# x2$samples$norm.factors <- 1
# x2$counts[,1] <- ceiling(x2$counts[,1]*0.05) # reduce first sample to 5%
# x2$counts[,2] <- x2$counts[,2]*5 # Inflate second sample by x5

### Graph
# # pdf("../Figures/Proof/theoretical_normalisation_effect.pdf")
# par(mfrow = c(1,2))
# lcpm <- cpm(x2, log = T)
# boxplot(lcpm, las = 2, col = col, main = "")
# title(main = "A. Example: Unnormalised data", ylab = "Log-cpm")
# x2 <- calcNormFactors(x2)
# x2$samples$norm.factors
# lcpm <- cpm(x2, log = T)
# boxplot(lcpm, las = 2, col = col, main = "")
# title(main = "B. Example: Normalised data", ylab = "Log-cpm")
# # dev.off()

# Unsupervised clustering of samples & runs
# lcpm <- cpm(x, log = T)
# par(mfrow = c(1,2))
# col.group <- as.factor(group)
# levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
# col.group <- as.character(col.group)
# col.lane <- as.factor(lane)
# levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
# col.lane <- as.character(col.lane)
# pdf("../Figures/Proof/Groups_MDS.pdf")
# par(mfrow = c(1, 1))
# plotMDS(lcpm, labels = colnames(x), col = col.group)
# title(main = "A. Sample groups")
# dev.off()

# plotMDS(lcpm, labels = lane, col = col.lane, dim = c(1,2))
# title(main = "B. Sequencing lanes")

### Double check...
lcpm <- cpm(x, log = T)

## Online launch of this.
# biocLite("Glimma", dependencies = T)
# library(Glimma)
# glMDSPlot(lcpm, labels = paste(colnames(x), sep = "_"),
#           groups = x$samples[,c(2)], launch = T, top = 13564)
# col.cell <- c("#999999","#56B4E9","#E69F00","#009E73","#CC79A7")[x$samples$group]
# data.frame(x$samples$group, col.cell)
# pdf("../Figures/Proof/PCA_of_all_genes.pdf")
col.cell <- c("#999999","#56B4E9","#E69F00","#009E73","#CC79A7")[x$samples$group]
?plotMDS

plotMDS(lcpm, pch = 16, cex = 2, col = col.cell, top = 12616, labels = colnames(x))
legend("top",
       fill = c("#999999", "#56B4E9",
                "#E69F00", "#009E73",
                "#CC79A7"),
       legend = levels(x$samples$group))
# Add a title
# title(main = "A. PCA of log(cpm) for all genes")
# dev.off()

# Get the gene names for the top 100 most variable genes
## Estimate variance for all genes
var_genes <- apply(lcpm, 1, var)

## Remove TRAV and TRBV
att <- as.data.frame(var_genes)
att1 <- rownames_to_column(att, var = "ENTREZID")
att2 <- merge(att1, x$genes, by = "ENTREZID")
att3 <- att2[!grepl("^TRAV", att2$SYMBOL),]
att4 <- att3[!grepl("^TRBV", att3$SYMBOL),]
var_genes2 <- extract2(att4, 'var_genes') %>% set_names(att4$ENTREZID)

select_var <- names(sort(var_genes2, decreasing = T))[1:150]

# Subset logcounts matrix
highly_variable_lcpm <- lcpm[select_var,]
# dim(highly_variable_lcpm)
# head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycol <- colorpanel(1000,"blue","white","red")

## Get the SYMBOLS instead of ENTREZ ID
this <- rownames_to_column(as.data.frame(highly_variable_lcpm), var = "ENTREZID")
this1 <- merge(this,x$genes, by = "ENTREZID")
highly_variable_lcpm_sym <- within(this1, rm(TXCHROM, ENTREZID))
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

# Plot the heatmap
# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/Heatmap_top150_variable_genes.png",
#     width = 300, height = 350, units = "mm", res = 300)
# heatmap.2(as.matrix(highly_variable_lcpm_sym[, names(highly_variable_lcpm_sym) != "SYMBOL"]),
#           col = mycol,
#           trace = "none",
#           density.info = "none",
#           main = "B. Heatmap of the Top 150 most Variable Genes",
#           cex.main = 1.5,
#           #ColSideColors = col.cell, 
#           scale = "row",
#           margin = c(10,5), lhei = c(2,10),
#           labCol = colnames(highly_variable_lcpm_sym),
#           labRow = highly_variable_lcpm_sym$SYMBOL,
#           hclustfun = hclustAvg)
# par(xpd = T)
# legend(x = 0.87, y = 1.065,
#        fill = c("#999999", "#56B4E9",
#                 "#E69F00", "#009E73",
#                 "#CC79A7"),
#        legend = levels(x$samples$group))
# dev.off()


# EntrezID as row
# heatmap.2(highly_variable_lcpm, 
#           col = mycol,
#           trace = "none", 
#           main = "Top 500 most variable genes across samples",
#           ColSideColors = col.cell,scale = "row",
#           margin = c(10,10), lhei = c(2,10),
#           labCol = colnames(highly_variable_lcpm),
#           labRow = rownames(highly_variable_lcpm))

# Differential gene expression
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  CD27LOvsCD27HI = VD1.CD27LO - VD1.CD27HI, 
  CD27HIvsEMRA = VD1.CD27HI - CD8.EMRA, 
  CD27HIvsNaive = VD1.CD27HI - CD8.Naive, 
  CD27HIvsVD2 = VD1.CD27HI - VD2,
  CD27LOvsEMRA = VD1.CD27LO - CD8.EMRA,
  CD27LOvsNaive = VD1.CD27LO - CD8.Naive,
  CD27LOvsVD2 = VD1.CD27LO - VD2,
  EMRAvsNaive = CD8.EMRA - CD8.Naive,
  EMRAvsVD2 = CD8.EMRA - VD2,
  NaivevsVD2 = CD8.Naive - VD2,
  levels = colnames(design))

## Removing heteroscedascity
# biocLite("limma", dependencies = T)
library(limma)
# pdf("../Figures/Proof/Mean-variance-tred.pdf")
# par(mfrow = c(1,2))
v <- voom(x, design, plot = F)
# save.image(file = "half.RData")
# load("Bulk/Counts/half.RData")

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = T)
# plotSA(efit, main = "Final model: Mean−variance trend")
# dev.off()

## Number of differentially expressed genes
summary(decideTests(efit))

### More stringently selected DE genes
tfit <- treat(vfit, lfc = 1)
summary(decideTests(tfit))

dt <- decideTests(efit)

head(dt)
# summary(dt)

de.common <- which(dt[,1]!=0 & dt[,8]!=0)
# length(de.common)

#### Write the common genes
Common_genes <- efit$genes$ENTREZ[de.common]
# Entrez <- tfit$genes$ENTREZ[de.common]
# Symbol <- tfit$genes$SYMBOL[de.common]
# Common_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(Common_DE_genes)

VD1  <- which(dt[,1]!=0 & dt[,8]==0)
# length(VD1)
de.VD1 <- efit$genes$ENTREZID[VD1]
# Entrez <- tfit$genes$ENTREZ[de.VD1]
# Symbol <- tfit$genes$SYMBOL[de.VD1]
# VD1_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(VD1_DE_genes)

CD8 <- which(dt[,1]==0 & dt[,8]!=0)
de.CD8 <- efit$genes$ENTREZID[CD8]
# Entrez <- tfit$genes$ENTREZ[de.CD8]
# Symbol <- tfit$genes$SYMBOL[de.CD8]
# CD8_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(CD8_DE_genes)
# http://eulerr.co/


## Up, Down and Both
# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/DE_genes_EffectorvsNaive.png",
#     width = 450, height = 300, units = "mm", res = 300)
# layout(matrix(c(1,1,2,3), 2, 2, byrow = F))
# vennDiagram(dt[, c(1,8)], include = "both", mar = c(0,0,0,0))
# title(main = "Differentially expressed Genes", line = -9)
# title(main = "C. Differential expression between effector and naive subsets", line = -2, cex.main = 1.5)
# vennDiagram(dt[, c(1,8)], include = "up", mar = c(0,0,0,0))
# title(main = "Upregulated", line = -2)
# vennDiagram(dt[, c(1,8)],  include = "down", mar = c(0,0,0,0))
# title(main = "Downregulated", line = -2)
# dev.off()

# Heatmap of shared differentially expressed genes 
## Remove VD2
v1 <- v
check <- as.data.frame(v$E)
v1$targets <- droplevels(subset(v1$targets, group != "VD2"))
v1$E <- as.data.frame(v1$E)
v1$E <- v1$E %>% dplyr::select(-contains("VD2"))
v1$E <- as.matrix(v1$E)
j <- !grepl("VD2", group)
group1 <- group[j]
i <- which(v1$genes$ENTREZID %in% Common_genes)


mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v1$E[i,], scale = "row",
#           labRow = v1$genes$SYMBOL[i], labCol = colnames(v1), #can also do labCol = groups
#           col = mycol, trace = "none", density.info = "none",
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "Shared DE genes (Naive v Effectors)",
#           hclustfun = hclustAvg)
# dev.off()

# PCA of all shared differentially expressed genes
## Change all 'v' to 'v1' to remove VD2 from PCA
library(ggbiplot)
# i <- which(v1$genes$ENTREZID %in% Common_genes)
# this <- as.data.frame(v1$E[i,])
# that <- rownames_to_column(this, var = "Gene")
# head(that)
# 
# that1 <- that %>% gather(-contains("Gene"), key = "Sample", value = "value")
# that2 <- spread(that1, key = Gene, value = value)
# 
# # Label groups
# subgroup <- ifelse(grepl("VD1.CD27LO", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", that2$Sample), "VD1.CD27LO",
#                    ifelse(grepl("CD8.EMRA", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", that2$Sample), "CD8.EMRA",
#                           ifelse(grepl("CD8.Naive", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", that2$Sample), "CD8.Naive",
#                                  ifelse(grepl("VD1.CD27HI", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", that2$Sample), "VD1.CD27HI",
#                                         ifelse(grepl("VD2", that2$Sample), "VD2", "none")))))
# that2$subgroup <- as.factor(subgroup)

# Remove identifiers
# that3 <- data.frame(that2[, names(that2) != "Sample"], row.names = that2[, names(that2) == "Sample"])
# colnames(that3) <- gsub("^X", "",  colnames(that3))
# 
# prin_comp <- prcomp(that3[, names(that3) != "subgroup"], scale. = T)
# subgroup <- that3[, "subgroup"]
# 
# dev.off()
# g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
#               groups = subgroup, ellipse = F,
#               circle = T,
#               var.axes = F
# )
# g <- g + scale_color_manual(values = cbcols)
# g <- g + theme_bw()
# g <- g + theme(legend.direction = 'horizontal', 
#                legend.position = 'top')
# 
# g <- g + ggtitle("PCA of shared DE genes (Naive vs Effector)")
# 
# ggsave("PCA of shared DE genes (Naive vs Effector).png" , plot = g, device = "png",
#        path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures",
#        height = 6, width = 6, units = 'in', dpi = 600)

# # How many PC's are important?
# library(factoextra)
# fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue")
# 
# # Find the Eigenvalues
# eig <- (prin_comp$sdev)^2
# 
# ## Variances in percentage
# variance <- eig*100/sum(eig)
# 
# ## Cumulative variances
# cumvar <- cumsum(variance)
# 
# ## Store the variances as a dataframe
# ### Write as a DF
# eigenvalues_RNAseq <- data.frame(eig = eig, variance = variance,
#                                cumvariance = cumvar)
# # writeCsvO(eigenvalues_RNAseq)
# 
# # Store the variances
# var <- get_pca_var(prin_comp)
# 
# ## Find the coordinates of variables
# ### var$coord[, 1:5]
# 
# ## Find the correlation between variables and principal components
# loadings <- prin_comp$rotation
# 
# # ### Find orientation of loadings
# loads <- as.data.frame(loadings)
# loads[order(loads$PC1, decreasing = T)[1:10],]
# loads1 <- tibble:: rownames_to_column(loads, "Parameter")
# loads_contrib <- droplevels(subset(loads1, PC1 < -0.26))
# loads_contrib$Parameter <- as.factor(loads_contrib$Parameter)
# 
# sdev <- prin_comp$sdev
# 
# ### Find the correlataions
# var.coord <- t(apply(loadings, 1, var_cor_func, sdev))
# 
# ## Calculate the Cos2 (square of the coordinates)
# var.cos2 <- var.coord^2
# 
# ## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
# comp.cos2 <- apply(var.cos2, 2, sum)
# var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
# contrib.var <- as.data.frame(var.contrib)
# 
# ## Find the most contributing variable
# contrib.var[order(contrib.var$PC1, decreasing = T)[1:10],]
# 
# # Use the analysis to compare the values for the groups across PC
# # Store the contribution to each PC
# PC_TCGA <- as.data.frame(prin_comp$x)
# PC_TCGA1 <- tibble:: rownames_to_column(PC_TCGA, "Sample")
# PC_TCGA1$Sample <- as.factor(PC_TCGA1$Sample)
# 
# ## Gather the PCs 
# PC_TCGA2 <- PC_TCGA1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
# PC_TCGA2$Component <- as.factor(PC_TCGA2$Component)
# 
# # Plot
# ## Test for normality
# subgroup <- ifelse(grepl("VD1.CD27LO", PC_TCGA2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", PC_TCGA2$Sample), "VD1.CD27LO",
#                    ifelse(grepl("CD8.EMRA", PC_TCGA2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", PC_TCGA2$Sample), "CD8.EMRA",
#                           ifelse(grepl("CD8.Naive", PC_TCGA2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", PC_TCGA2$Sample), "CD8.Naive",
#                                  ifelse(grepl("VD1.CD27HI", PC_TCGA2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", PC_TCGA2$Sample), "VD1.CD27HI",
#                                         ifelse(grepl("VD2", PC_TCGA2$Sample), "VD2", "none")))))
# PC_TCGA2$subgroup <- as.factor(subgroup)
# Normality <- PC_TCGA2
# Normality$uniq <- as.factor(paste(Normality$subgroup, Normality$Component, sep = ","))
# 
# ## Shapiro test
# normal_list <- list()
# c <- 1
# for(i in levels(Normality$uniq)){
#   name <- basename(i)
#   cat('Processing', i, '\n')
#   working <- droplevels(subset(Normality, uniq == i))
#   stat <- shapiro.test(working[,"ComponentScore"])
#   normal <- as.data.frame(stat$p.value)
#   normal_list[[i]] <- cbind(i, normal)
#   c <- c + 1
# }
# z <- do.call(rbind, normal_list)
# rownames(z) <- c()
# head(z)
# 
# ## QQ plots
# for(i in levels(Normality$uniq)){
#   name <- basename(i)
#   cat('Processing', i, '\n')
#   working <- droplevels(subset(Normality, uniq == i))
#   temp_plot <- ggplot(working) +
#     stat_qq(aes(sample = ComponentScore))
#   filen <- paste0(i, ".png")
#   ggsave(filen, plot = temp_plot, device = "tiff",
#          path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/PCA/Normality_PC",
#          height=5, width=5, units='in', dpi=600)
# }
# 
# # Compare across components 
# for(i in levels(PC_TCGA2$Component)){
#   name <- basename(i)
#   cat('Processing', i, '\n')
#   Chosen <- droplevels(subset(PC_TCGA2, Component == i))
#   Chosen$Rank <- rank(Chosen$ComponentScore)
#   temp_plot <- ggComponent(Chosen)
#   filen <- paste0(i, ".png")
#   ggsave(filen, plot = temp_plot, device = "png",
#          path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/PCA/Compare_PC",
#          height=6, width=7, units='in', dpi=600)
# }
# 
# # Combine stats into a list
# stat_list <- list()
# c <- 1
# for(i in levels(PC_TCGA2$Component)) {
#   name <- basename(i)
#   cat('Processing', i, '\n')
#   workingon <- droplevels(subset(PC_TCGA2, Component == i))
#   workingon$Rank <- rank(workingon$ComponentScore)
#   # assign your ggplot call to the i'th position in the list
#   yz <- compare_means(Rank ~ subgroup, data = workingon, method = "wilcox.test")
#   y <- as.data.frame(yz)
#   stat_list[[i]]  <- cbind(i, y)
#   c <- c + 1}
# 
# # Bind and remove row names
# z <- do.call(rbind, stat_list)
# head(z)
# rownames(z) <- c()
# RNAseq_PCAStats <- z

# Write out the statistics
# writeCsvO(RNAseq_PCAStats)

# write.fit(tfit, dt, file="results.txt")


# View comparisons
CD27LO.vs.CD27HI <- topTreat(efit, coef = 1, n = Inf)
CD27HI.vs.EMRA <- topTreat(efit, coef = 2, n = Inf)
CD27HI.vs.Naive <- topTreat(efit, coef = 3, n = Inf)
CD27HI.vs.VD2 <- topTreat(efit, coef = 4, n = Inf)

CD27LO.vs.EMRA <- topTreat(efit, coef = 5, n = Inf)
CD27LO.vs.Naive <- topTreat(efit, coef = 6, n = Inf)
CD27LO.vs.VD2 <- topTreat(efit, coef = 7, n = Inf)

EMRA.vs.Naive <- topTreat(efit, coef = 8, n = Inf)
EMRA.vs.VD2 <- topTreat(efit, coef = 9, n = Inf)

Naive.vs.VD2 <- topTreat(efit, coef = 10, n = Inf)

# str(CD27LO.vs.CD27HI)
sig <- droplevels(subset(CD27LO.vs.CD27HI, adj.P.Val <= 0.05))
nrow(sig)


sig <- droplevels(subset(CD27LO.vs.EMRA, adj.P.Val <= 0.05))
nrow(sig)


CD27HI.vs.Naive[grepl("TRAV", CD27HI.vs.Naive$SYMBOL), ]
CD27LO.vs.EMRA[grepl("TRAV", CD27LO.vs.EMRA$SYMBOL), ]

CD27HI.vs.Naive[grepl("TRDV", CD27HI.vs.Naive$SYMBOL), ]



save.image("../Final_edgeR.RData")

# These are the outputs of EdgeR!!
# head(CD27LO.vs.EMRA)
# head(CD27HI.vs.VD2)
# head(CD27HI.vs.EMRA)
# head(CD27HI.vs.Naive)

# plotMD(tfit, column = 1, status = dt[,1], main = colnames(efit)[1], 
#        xlim = c(-8,13))
# 
# glMDPlot(tfit, coef = 1, status = dt, main = colnames(efit)[1],
#          side.main = "ENTREZID", counts = x$counts, groups = group, launch = T)

# Heatmaps
## Change the top genes, for various heatmaps looking at different sets
library(gplots)

# # CD27HI versus CD27LO
# CD27LO.vs.CD27HI1 <- CD27LO.vs.CD27HI[!grepl("^TRAV", CD27LO.vs.CD27HI$SYMBOL),]
# CD27LO.vs.CD27HI2 <- CD27LO.vs.CD27HI1[!grepl("^TRBV", CD27LO.vs.CD27HI1$SYMBOL),]
# CD27LO.vs.CD27HI2 <- Take_Sigs(CD27LO.vs.CD27HI2)
# 
# CD27LO.vs.CD27HI.topgenes <- CD27LO.vs.CD27HI2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.CD27HI.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs CD27HI")

# gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.CD27HI2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "CD27HI", "CD27LO") # THIS WILL CHANGE ONCE YOU SET THE LFC = 2
# that <- droplevels(subset(this, Group == "CD27HI"))
# 
# CD27LO_CD27HI_DEgenes <- this
# writeCsvO(CD27LO_CD27HI_DEgenes)

# # CD27HI versus EMRA
# CD27HI.vs.EMRA1 <- CD27HI.vs.EMRA[!grepl("^TRAV", CD27HI.vs.EMRA$SYMBOL),]
# CD27HI.vs.EMRA2 <- CD27HI.vs.EMRA1[!grepl("^TRBV", CD27HI.vs.EMRA1$SYMBOL),]
# CD27HI.vs.EMRA2 <- Take_Sigs(CD27HI.vs.EMRA2)
# 
# CD27HI.vs.EMRA.topgenes <- CD27HI.vs.EMRA2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.EMRA.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27HI vs EMRA")
## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27HI")
# 
# CD27HI_EMRA_DEgenes <- this
# writeCsvO(CD27HI_EMRA_DEgenes)

# # CD27HI versus Naive
# CD27HI.vs.Naive1 <- CD27HI.vs.Naive[!grepl("^TRAV", CD27HI.vs.Naive$SYMBOL),]
# CD27HI.vs.Naive2 <- CD27HI.vs.Naive1[!grepl("^TRBV", CD27HI.vs.Naive1$SYMBOL),]
# CD27HI.vs.Naive2 <- Take_Sigs(CD27HI.vs.Naive2)
# 
# CD27HI.vs.Naive.topgenes <- CD27HI.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.Naive.topgenes)
# 
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27HI vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27HI")
# CD27HI_Naive_DEgenes <- this
# writeCsvO(CD27HI_Naive_DEgenes)


# # CD27HI versus VD2
# CD27HI.vs.VD21 <- CD27HI.vs.VD2[!grepl("^TRAV", CD27HI.vs.VD2$SYMBOL),]
# CD27HI.vs.VD22 <- CD27HI.vs.VD21[!grepl("^TRBV", CD27HI.vs.VD21$SYMBOL),]
# CD27HI.vs.VD22 <- Take_Sigs(CD27HI.vs.VD22)
# 
# CD27HI.vs.VD2.topgenes <- CD27HI.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27HI vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27HI")
# CD27HI_VD2_DEgenes <- this
# writeCsvO(CD27HI_VD2_DEgenes)


# # CD27LO versus VD2
# CD27LO.vs.VD21 <- CD27LO.vs.VD2[!grepl("^TRAV", CD27LO.vs.VD2$SYMBOL),]
# CD27LO.vs.VD22 <- CD27LO.vs.VD21[!grepl("^TRBV", CD27LO.vs.VD21$SYMBOL),]
# CD27LO.vs.VD22 <- Take_Sigs(CD27LO.vs.VD22)
# 
# CD27LO.vs.VD2.topgenes <- CD27LO.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27LO")
# CD27LO_VD2_DEgenes <- this
# writeCsvO(CD27LO_VD2_DEgenes)


# # CD27LO versus EMRA
# CD27LO.vs.EMRA1 <- CD27LO.vs.EMRA[!grepl("^TRAV", CD27LO.vs.EMRA$SYMBOL),]
# CD27LO.vs.EMRA2 <- CD27LO.vs.EMRA1[!grepl("^TRBV", CD27LO.vs.EMRA1$SYMBOL),]
# CD27LO.vs.EMRA2 <- Take_Sigs(CD27LO.vs.EMRA2)
# 
# CD27LO.vs.EMRA.topgenes <- CD27LO.vs.EMRA2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.EMRA.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs EMRA")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27LO")
# CD27LO_EMRA_DEgenes <- this
# writeCsvO(CD27LO_EMRA_DEgenes)


# # CD27LO versus Naive
# CD27LO.vs.Naive1 <- CD27LO.vs.Naive[!grepl("^TRAV", CD27LO.vs.Naive$SYMBOL),]
# CD27LO.vs.Naive2 <- CD27LO.vs.Naive1[!grepl("^TRBV", CD27LO.vs.Naive1$SYMBOL),]
# CD27LO.vs.Naive2 <- Take_Sigs(CD27LO.vs.Naive2)
# 
# CD27LO.vs.Naive.topgenes <- CD27LO.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.Naive.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27LO vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27LO")
# CD27LO_Naive_DEgenes <- this
# writeCsvO(CD27LO_Naive_DEgenes)


# # EMRA versus Naive
# EMRA.vs.Naive1 <- EMRA.vs.Naive[!grepl("^TRAV", EMRA.vs.Naive$SYMBOL),]
# EMRA.vs.Naive2 <- EMRA.vs.Naive1[!grepl("^TRBV", EMRA.vs.Naive1$SYMBOL),]
# EMRA.vs.Naive2 <- Take_Sigs(EMRA.vs.Naive2)
# 
# EMRA.vs.Naive.topgenes <- EMRA.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% EMRA.vs.Naive.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "EMRA vs Naive")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "EMRA")
# EMRA_Naive_DEgenes <- this
# writeCsvO(EMRA_Naive_DEgenes)


# # EMRA versus VD2
# EMRA.vs.VD21 <- EMRA.vs.VD2[!grepl("^TRAV", EMRA.vs.VD2$SYMBOL),]
# EMRA.vs.VD22 <- EMRA.vs.VD21[!grepl("^TRBV", EMRA.vs.VD21$SYMBOL),]
# EMRA.vs.VD22 <- Take_Sigs(EMRA.vs.VD22)
# 
# EMRA.vs.VD2.topgenes <- EMRA.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% EMRA.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "EMRA vs VD2")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "EMRA")
# EMRA_VD2_DEgenes <- this
# writeCsvO(EMRA_VD2_DEgenes)


# # Naive versus VD2
# Naive.vs.VD21 <- Naive.vs.VD2[!grepl("^TRAV", Naive.vs.VD2$SYMBOL),]
# Naive.vs.VD22 <- Naive.vs.VD21[!grepl("^TRBV", Naive.vs.VD21$SYMBOL),]
# Naive.vs.VD22 <- Take_Sigs(Naive.vs.VD22)
# 
# Naive.vs.VD2.topgenes <- Naive.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% Naive.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "Naive vs VD2")

## gaining a dataframe with the differentially expressed
# this <- Naive.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "Naive")
# Naive_VD2_DEgenes <- this
# writeCsvO(Naive_VD2_DEgenes)


# # Shared top 100 DE genes
# k <- which(CD27LO.vs.CD27HI.topgenes %in% EMRA.vs.Naive.topgenes)
# j <- CD27LO.vs.CD27HI.topgenes[k]
# i <- which(v1$genes$ENTREZID %in% j)
# mycol <- colorpanel(1000,"blue","white","red")
# 
# col.cell1 <- c("#999999","#56B4E9","#E69F00","#009E73")[v1$targets$group]
# data.frame(v1$targets$group, col.cell1)
# 
# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/Shared_top100_DEgenes.png",
#     width = 300, height = 300, units = "mm", res = 300)
# heatmap.2(v1$E[i,], scale = "row",
#           labRow = v1$genes$SYMBOL[i], labCol = group1, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), ColSideColors = col.cell1,
#           main = "Shared top 100 DE genes (Naive v Effectors)")
# par(xpd = T)
# legend(x = 0.87, y = 1.05, 
#        fill = c("#999999", "#56B4E9",
#                 "#E69F00", "#009E73"),
#        legend = levels(v1$targets$group))
# dev.off()


# # Check other way for %in% (same)
# k <- which(EMRA.vs.Naive.topgenes  %in% CD27LO.vs.CD27HI.topgenes)
# j <- EMRA.vs.Naive.topgenes[k]
# i <- which(v1$genes$ENTREZID %in% j)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v1$E[i,], scale = "row",
#           labRow = v1$genes$SYMBOL[i], labCol = group1, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "Shared top 100 DE genes (Naive v Effectors)")



# Find top 10 genes in all combinations
# colnames(contr.matrix)
# 
# this <- as.vector(rbind(CD27LO.vs.CD27HI2$ENTREZID[1:25], 
#                         CD27HI.vs.EMRA2$ENTREZID[1:25], 
#                         CD27HI.vs.Naive2$ENTREZID[1:25],
#                         CD27HI.vs.Naive2$ENTREZID[1:25],
#                         CD27HI.vs.VD22$ENTREZID[1:25],
#                         CD27LO.vs.EMRA2$ENTREZID[1:25],
#                         CD27LO.vs.Naive2$ENTREZID[1:25],
#                         CD27LO.vs.VD22$ENTREZID[1:25],
#                         EMRA.vs.Naive2$ENTREZID[1:25],
#                         EMRA.vs.VD22$ENTREZID[1:25],
#                         Naive.vs.VD22$ENTREZID[1:25]))
# this <- this[!duplicated(this)]
# 
# i <- which(v$genes$ENTREZID %in% this)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "Top 20 DE genes across all populations")
# length(this)

#### Write the common genes








# Geneset Enrichment Analysis
## GO_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/GO_genesets.rdata")
GO_sizes<- as.data.frame(lengths(Hs.c5)) %>% rownames_to_column(., var = "GO_pathway")
colnames(GO_sizes) <- c("GO_Pathway", "Size")
remove_these <- droplevels(subset(GO_sizes, Size >= 300 & Size <= 10))$GO_Pathway
isNameInIndex <- names(Hs.c5) %in% remove_these
GO_terms <- Hs.c5[!isNameInIndex]

idx <- ids2indices(GO_terms, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.05))


att <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "GenesetName")
head(att)

these <- att$GenesetName

isNameInIndex <- names(GO_terms) %in% these
these1 <- GO_terms[isNameInIndex]

head(these1)

GO_Genesets_Sig_GD <- do.call(rbind, these1)

GO_Genesets_Sig_GD$id <- rep(names(these1), sapply(these1, nrow))


write
barcodeplot(efit$t[, 1], index = idx$GO_CELL_KILLING, main = "CD27LO vs CD27HI")

# Immunological Signatures
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Immunological_Signatures.rdata")
idx <- ids2indices(Hs.c7, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 10)

barcodeplot(efit$t[, 1], index = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP , 
index2 = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN, main = "NaiveVsCD8")

# Hallmark_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Hallmark_genesets.rdata")
idx <- ids2indices(Hs.H, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, PValue < 0.05))

# KEGG
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/kegg_human.rdata")
idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)

getwd()
write.table(cam.CD27LO.vs.CD27HI, file = "../Output/Cyto/CD27LO.CD27HI_Genesets.txt")

idx <- ids2indices(kegg_human, id = rownames(v))
cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)


idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.EMRA <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsEMRA"])
head(cam.CD27LO.vs.EMRA, 5)


idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27HI.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27HIvsNaive"])
head(cam.CD27HI.vs.Naive, 5)

View(contr.matrix)
# Searching
## CD8 diffs
this <- rownames_to_column(cam.EMRA.vs.Naive, var = "KEGG")
droplevels(subset(this, KEGG == "hsa00190 Oxidative phosphorylation"))

that <- droplevels(subset(this, FDR < 0.01))

that1 <- that$KEGG

this <- rownames_to_column(cam.CD27LO.vs.CD27HI, var = "KEGG")
that <- droplevels(subset(this, FDR < 0.01))
length(that$KEGG)

 
