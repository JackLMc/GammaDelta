# A script to investigate the single-cell RNA sequencing data from Mike Stubbington on GammaDelta T-cells
# Original Author: Xi Chen (pre Seurat Version 2.3.4)
# Co-Author: Jack McMurray (Updated to Seurat 2.3.4, and then to version 4.0.0)

###############################################################################
#### This is the final set of transcriptomic analyses of the gamma-delta cells. 
#### This uses only the cells that passed QC based on read numbers etc. 
#### It's also based on the initial exploratory analyses done elsewhere.
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
load("RData/SeuratAnalyses.RData")

# #### Read in the table of TPM expression values for the cells that passed QC - 
# ##### This was performed outside of R
# tpms <- read.table("Data/clean_quant_table.csv", sep = ",", header = T, row.names = 1)
# 
# #### Rename the gene names from ENSEMBL identifiers to gene symbols to make it easier to interpret results
# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- row.names(tpms)
# G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"),
#                 values = genes, mart = mart)
# 
# backup_g_list <- G_list
# G_list <- backup_g_list
# row.names(G_list) <- G_list$ensembl_gene_id
# G_list <- G_list[, 2, drop = F] 
# 
# 
# dim(tpms)
# tpms <- tpms[row.names(G_list), ]
# dim(tpms)
# 
# #### We now have a TPM table with gene names instead of ENSEMBL IDs. 
# #### We've lost 1434 rows because they're strange entries which aren't really genes.
# #### This is nearly double what was lost in original analysis (620)
# #### We can safely ignore these but it's better to rescale the columns so they sum to 1 million.
# rescale <- function(x) {
#     s = sum(x)
#     (x/s) * 1e6}
# 
# tpms <- apply(tpms, 2, rescale)
# 
# #### There are cases where more than one ENSEMBL ID resolves to the same gene name
# #### We want to add the TPMs for each duplicated name.
# library(plyr)
# 
# tpms_with_names <- merge(tpms, G_list, by = 0)
# tpms_with_names <- subset(tpms_with_names, select = -c(Row.names))
# tpms_with_names <- plyr:: ddply(tpms_with_names, "external_gene_name", numcolwise(sum))
# 
# rownames(tpms_with_names) <- tpms_with_names$external_gene_name
# tpms <- subset(tpms_with_names, select = -c(external_gene_name))
# 
# ## Setting up Seurat object ##
# #### Normalisation is for Counts, so avoid doing that and pass in log transformed TPMs
# tpms <- log(tpms + 1)
# 
# #### Check how many genes have at least one transcript in each cell
# at_least_one <- apply(tpms, 2, function(x) sum(x>0))
# hist(at_least_one, breaks = 100,
#      main = "Distribution of detected genes",
#      xlab = "Genes with at least one tag")
# 
# hist(colSums(tpms),
#      breaks = 100, main = "Expression sum per cell",
#      xlab = "Sum expression")
# 
# #### Manually check the number of genes detected in three or more cells
# #### A lot of genes are not detected in 3 or more cells
# tmp <- apply(tpms, 1, function(x) sum(x>0))
# table(tmp>=3)
# #### all cells have at least 2000 detected genes (Min)
# 
# keep <- tmp>=3
# tmp <- tpms[keep, ]
# at_least_one <- apply(tmp, 2, function(x) sum(x>0))
# summary(at_least_one)
# 
# #### Create the object with new Seurat workflow
# #### Using the parameters that Mike used in his script (not informed from above)
# #### Can't be done using `Setup` anymore, function is deprecated
# gd.data <- CreateSeuratObject(counts = tpms,
#                               min.cells = 3,
#                               min.genes = 1000,
#                               is.expr = 1,
#                               names.field = 2,
#                               names.delim = "_",
#                               project = "gammadelta")
# 
# ## Check Mitochondrial genes ##
# #### A good QC metric is the percentage of mitochrondrial RNA
# #### Cells that were lysed, will lose cytoplasmic RNA, but mito will be kept and sequenced
# gd.data[["percent.mt"]] <- PercentageFeatureSet(gd.data, pattern = "^MT-")
# 
# VlnPlot(gd.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 
# gd.data <- subset(gd.data, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5) # Choose these based on previous plots
# 
# #### Not a problem, all cells have below 5% mitochondrial RNA
# 
# gd.data <- FindVariableFeatures(gd.data, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(gd.data), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(gd.data)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
# plot1 + plot2
# 
# 
# #### Although, initially, I thought that there weren't batch effects in the data.
# #### It turns out that there are and that Donor 6 sits apart from the others.
# #### This makes sense since that plate was done at a very different time.
# #### So, I want to remove the batch effects (this is nice and easy because the plates aren't confounded).
# #### Exploratory analyses found that the number of genes detected influences the results and that the biology appears easier to interpret if we remove this.
# all.genes <- rownames(gd.data)
# gd.data <- ScaleData(gd.data, features = all.genes, vars.to.regress = c("orig.ident", "nCount_RNA"))
# # Regressing out the orig.ident and the nCount_RNA as to remove any batch effects
# 
# # SKIP THE NORMALISATION STEP HERE. SINCE I PASSED IT LOG(TPM + 1)
# 
# #### Next we want to find a set of genes that are sufficiently variable
# #### (and expressed at sufficient levels) to be interesting.
# #### This is a *bit* arbritary but these cutoffs are fairly relaxed so I don't think we're missing much.
# #### Making them more relaxed tended to make the results look weird and not in an interesting way.

## Initial analyses ##
#### Run PCA
gd.data <- RunPCA(gd.data, features = VariableFeatures(object = gd.data))
print(gd.data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(gd.data, dims = 1:2, reduction = "pca")


#### There seems to be some outlier cells
# Find Outlier cells
DimPlot(gd.data, reduction = "pca")

#### There are outliers, which make it difficult to see the strucutre in the main population.
#### These are likely technical artefacts, so we'll remove them.
rotations <- gd.data@reductions$pca@cell.embeddings %>% as.data.frame()

outlier_cells <- rotations[(rotations$PC_1>(20) | rotations$PC_2>(20)), ]
tpms_without_outliers <- tpms[, !colnames(tpms) %in% rownames(outlier_cells)]

dim(tpms)
dim(tpms_without_outliers)

## 2 outlier cells have been removed

#### REPEAT WITH THE NONE OUTLIER CELLS ####
gd.data <- CreateSeuratObject(counts = tpms_without_outliers,
                              min.cells = 3,
                              min.genes = 1000,
                              is.expr = 1,
                              names.field = 2,
                              names.delim = "_",
                              project = "gammadelta")

## Check Mitochondrial genes ##
gd.data[["percent.mt"]] <- PercentageFeatureSet(gd.data, pattern = "^MT-")

VlnPlot(gd.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

gd.data <- subset(gd.data, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5) # Choose these based on previous plots
gd.data <- FindVariableFeatures(gd.data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gd.data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gd.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
plot1 + plot2

all.genes <- rownames(gd.data)
gd.data <- ScaleData(gd.data, features = all.genes, vars.to.regress = c("orig.ident", "nCount_RNA"))
# Regressing out the orig.ident and the nCount_RNA as to remove any batch effects

# SKIP THE NORMALISATION STEP HERE. SINCE I PASSED IT LOG(TPM + 1)

# Get the Stats
gd.data@meta.data
mean(gd.data@meta.data$percent.mt)

pdf("./Figures/CellReportsRewrite/SuppStatistics.pdf")
VlnPlot(gd.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
dev.off()


# Initial analyses again
gd.data <- RunPCA(gd.data, features = VariableFeatures(object = gd.data))
print(gd.data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gd.data, dims = 1:2, reduction = "pca")

### Produce some heatmaps about the data
DimHeatmap(gd.data, dims = 1, cells = 448, balanced = T)

#### Find Dimensionality of the data
gd.data <- JackStraw(gd.data, num.replicate = 100)
gd.data <- ScoreJackStraw(gd.data, dims = 1:20)


pdf("./Figures/CellReportsRewrite/JackStrawPlot.pdf") # Figure 1A
JackStrawPlot(gd.data, dims = 1:10)
dev.off()

# PC1 and PC2
### Calculate variance explained
pca <- gd.data@reductions$pca

mat <- Seurat::GetAssayData(gd.data, assay = "RNA", slot = "scale.data")
pca <- gd.data[["pca"]]

# Get the total variance:
eigValues <- (gd.data[["pca"]]@stdev)^2  ## EigenValues
varExplained <- eigValues / gd.data[["pca"]]@misc$total.variance ### VERY VERY LOW??

(gd.data[["pca"]]@stdev)^2/sum((gd.data[["pca"]]@stdev)^2) * 100



DimPlot(gd.data, reduction = "pca", cols = c("#999999", "#009E73", "#000000"), dims = c(1, 2))

gd.data <- FindNeighbors(gd.data, dims = 1:10)
gd.data <- FindClusters(gd.data, resolution = 0.5)

DimPlot(gd.data, reduction = "pca", cols = c("#999999", "#009E73"), dims = c(1, 2))


#### Run umap
# reticulate::py_install(packages = 'umap-learn')
gd.data <- RunUMAP(gd.data, dims = 1:10)
DimPlot(gd.data, reduction = "umap") # Doen't really improve it - try tSNE

#### Run tSNE
gd.data <- RunTSNE(gd.data, dims = 1:10)
DimPlot(gd.data, reduction = "tsne") # Doen't really improve it either - use PCA

#### Find all markers of Cluster 1
gd.data.markers <- FindAllMarkers(gd.data, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
gd.data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster0.markers <- FindMarkers(gd.data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = T)
cluster1.markers <- FindMarkers(gd.data, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = T)


## Making some plots
### Heatmap of genes contributing to PC1
# DimHeatmap(gd.data)

print(gd.data[["pca"]], dims = 1, nfeatures = 15)
PC1Features <- c("NKG7", "GZMH", "GZMB", "GNLY", "CST7", "CCL5", "GZMA", "CCL4", "TRBV11-3",
                 "FGFBP2", "CYTOR", "S100A4", "KLRF1", "CTSW", "CMC1", "RPL10A", "RPL30",
                 "SUSD3", "SNHG8", "RPS9", "RPL32", "RPS23", "EEF1G", "LDHB", "RPS13",
                 "SNHG29", "SELL", "NOSIP", "CCR7", "LTB")

# top10 <- gd.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(gd.data, features = PC1Features, group.colors = c("#999999", "#009E73"))


#### We calculated the PCAs using a subset of seemingly variable genes.
#### We now project the PCA analysis onto the entire dataset to give us weightings for all the genes.
#### This helps to ensure that we're not missing anything by filtering the genes
gd.data <- ProjectDim(gd.data)
DimHeatmap(gd.data, dims = 1) ### The PC1 heatmap is really nice. This is pretty much as we saw before
DimHeatmap(gd.data, dims = 2) ### PC2 isn't particularly interesting


## Interpreting the PCA ##
#### We want to check that the separation of cells on PC1 makes sense.
#### Lets see how the sort idenities (CD27+ or CD27-) map onto the plot
name_map <- read.csv("Data/name_map.txt", sep = "\t", row.names = 1)
rownames(name_map) = paste('X', rownames(name_map), sep = "")
name_map$sort_class <- ifelse((name_map$sort_class == "naive"), "VD1.CD27HI", "VD1.CD27LO")
name_map <- name_map[rownames(name_map) %in% rownames(gd.data[[]]), ]
name_map <- setNames(as.character(name_map$sort_class), rownames(name_map))
gd.data <- AddMetaData(object = gd.data, metadata = name_map, col.name = "sort_class")

DimPlot(gd.data, reduction = "pca", cols = c("#009E73", "#999999"), dims = c(1, 2), group.by = "sort_class")
DimPlot(gd.data, reduction = "pca", cols = c("#999999","#009E73"), dims = c(1, 2))








#### Do these split along the same lines as the sort classes
cluster.ids <- c("Cluster_1", "Cluster_2")
gd.data@meta.data$ident <- ifelse((gd.data@meta.data$seurat_clusters == "1"), "Cluster_1", "Cluster_2")

#### Plots to write out ####
data.plot <- gd.data@reductions$pca@cell.embeddings %>% as.data.frame()
cols_for_sort <- c("VD1.CD27LO" = "#999999", "VD1.CD27HI" = "#009E73")

data.plot_cluster <- merge(data.plot, gd.data@meta.data, by = "row.names", all.x = T)

pdf("./Figures/CellReportsRewrite/S1A_PCA_SortClass.pdf")
ggplot(data.plot_cluster, 
       mapping = aes(x = PC_1, y = PC_2, colour = sort_class)) +
    geom_point(shape = 17) +
    scale_color_manual(values = cols_for_sort) +
    theme(legend.title = element_blank()) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

### Clusters
cols_for_clusts <- c("Cluster_1" = "#009E73", "Cluster_2" = "#999999")

pdf("./Figures/CellReportsRewrite/1B_PCA_ClusteredClass.pdf")
ggplot(data.plot_cluster, 
       mapping = aes(x = PC_1, y = PC_2, colour = ident)) +
    geom_point(shape = 16) +
    scale_color_manual(values = cols_for_clusts) +
    theme(legend.title = element_blank()) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#### They are largely concordant!
#### There are some in the centre that aren't concordant, but this is either:
#### A limitation of the clustering
#### Representing the cells that were borderline of the FACS gates
#### The index sorting data would be really useful here to tease this apart - how concordant are they?
library(tidyverse)
clustered <- gd.data@meta.data %>% as.data.frame()

clustered$matches_ <- ifelse((clustered$ident == "Cluster_1"), "VD1.CD27HI", "VD1.CD27LO")
comparison <- clustered
total.cells <- nrow(comparison)

matching <- subset(comparison, sort_class == matches_)
matching.cells <- nrow(matching)

concordance <- matching.cells/total.cells
concordance # The sort and clustering phenotypes are 91% concordant
library(reshape2)
dcast(clustered, ident ~ sort_class)

continSort <- dcast(clustered, sort_class ~ ident) %>% column_to_rownames(., var  = "sort_class")
chisq_S <- chisq.test(continSort)


198 / (198 + 7) # Naive cells are 97% concordant with the cluster we thought
212 / (212 + 31) # Effector cells are 87% concordant with the cluster we thought

## Remake the Heatmap ##
#### Get scaled data, and the data to colour by
dat <- as.matrix(gd.data@assays$RNA@scale.data)

library(tidyverse)
meta <- as.data.frame(gd.data@meta.data) %>%
  rownames_to_column("Cells")
meta1 <- meta[, c("Cells", "ident")]
meta1$ident <- as.factor(meta1$ident)
col.cell1 <- c("#009E73", "#999999")[meta1$ident]
data.frame(meta1$ident, col.cell1)

#### Gain the topGenes
library(gplots)
HeatmapGenes <- c("NKG7", "GZMH", "GZMB", "GNLY", "CST7", "CCL5", "GZMA", "CCL4", "TRBV11-3",
                                  "FGFBP2", "CYTOR", "S100A4", "KLRF1", "CTSW", "CMC1", "RPL10A", "RPL30",
                                  "SUSD3", "SNHG8", "RPS9", "RPL32", "RPS23", "EEF1G", "LDHB", "RPS13",
                                  "SNHG29", "SELL", "NOSIP", "CCR7", "LTB")
genes.ordered <- rev(HeatmapGenes)

dat1 <- subset(dat, rownames(dat) %in% HeatmapGenes)
data.use <- as.data.frame(MinMax(data = dat1, min = -2.5, max = 2.5))

#### Better clustering
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

pdf("Figures/CellReportsRewrite/1C_PC1_Heat_plot.pdf", width = 3)
heatmap.2(as.matrix(data.use[genes.ordered,]),
          col = PurpleAndYellow(), trace = "none", cex.main = 1.5,
          Rowv = F, Colv = T, ColSideColors = col.cell1,
          scale = "none", margin = c(10,5), lhei = c(2,10),
          dendrogram = "none", hclustfun = hclustAvg, labCol = NA)
par(xpd = T)
legend(x = 0.87, y = 1.065,
       fill = c("#009E73", "#999999"),
       legend = levels(meta1$ident))
dev.off()

## Marker genes ##
#### Since PC1 appears to contain all the biology of interest, the best way to look for genes that are of interest
#### is to look at the PC1 weightings for the genes
#### Very negative weightings = EMRA-like cells, positive weightings = naive cells
#### How do some markers of interest look
pdf("Figures/CellReportsRewrite/1D_Mark_of_Interest.pdf", width = 18)
FeaturePlot(gd.data, c("NKG7","GZMB", "PRF1", "GNLY", "GZMA",
                       "CX3CR1", "LTB", "TCF7", "LEF1",
                       "IL7R", "NOSIP", "SELL"), reduction = "pca", ncol = 6)
dev.off()



#### As expected from the heatmaps, these look good

#### I looked at the genes with strong weightings in either direction in PC2 and they looked metabolic
weightings <- as.data.frame(gd.data@reductions$pca@feature.loadings)

PC1_orders <- weightings[order(weightings$PC_1), "PC_1", drop = F]
PC2_orders <- weightings[order(weightings$PC_2), "PC_2", drop = F]

write.table(PC1_orders, "Output/final_seurat_PC1_weightings.txt", sep = "\t", quote = F)
write.table(PC2_orders, "Output/final_seurat_PC2_weightings.txt", sep = "\t", quote = F)

sdev <- gd.data@reductions$pca@stdev
    
### Find the correlations
#### Functions
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}

var.coord <- t(apply(weightings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib) %>% rownames_to_column(var = "GeneID")

# colSums(contrib.var[, !'%in%'(colnames(contrib.var), "GeneID")])
# try <- contrib.var[contrib.var$GeneID %in% HeatmapGenes,]
# sum(try$PC_1)


# Produce a volcano plot of DGE
cluster.markers <- FindAllMarkers(gd.data, only.pos = F)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster.markers$Direction <- ifelse((cluster.markers$avg_log2FC > 0), "effector/Cluster2", "naive/Cluster1")

# head(cluster.markers, 30)
cluster.markers <- cluster.markers %>% rownames_to_column(., var  = "SYMBOL")
write.csv(cluster.markers, file = "./Output/SC_DE.csv") 


sigs <- droplevels(subset(cluster.markers, p_val_adj <= 0.05))


library(devtools)
# install_github("https://github.com/RGLab/MAST")

myData <- cluster.markers
myData$padjThresh <- as.factor(myData$p_val_adj < 0.05)
myData$Labels <- myData$SYMBOL
labelled_genes <- c("CX3CR1", "GZMA", "PRF1", "SELL", "IL7R", "CCR7", "TCF7", "LEF1", "EOMES", "TBX21")

myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""

myData[myData$SYMBOL %in% labelled_genes, ]

pdf(file = "./Figures/CellReportsRewrite/1E_Volcano.pdf", width = 6, height = 6)
library(ggrepel)
ggplot(data = myData, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(alpha = 0.2, size = 1) +
  theme(legend.position = "none", text = element_text(size = 10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ## Colours for significance
  geom_point(data = subset(myData, padjThresh==T & avg_log2FC<(-0.5)), aes(avg_log2FC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "blue4") +
  geom_point(data = subset(myData, padjThresh==T & avg_log2FC>0.5), aes(avg_log2FC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "orangered3") +
  ## Label these
  geom_label_repel(mapping = aes(label = Labels), min.segment.length = unit(0, "lines"), box.padding = 0.35, point.padding = 0.2, max.overlaps = 100) +
  # geom_text_repel(data=subset(myData,padjThresh==TRUE & avg_log2FC>0.5), aes(avg_log2FC,-log10(p_val),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
  labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
dev.off()

# #### Write out some tables that I want in other analyses in Python - Xi Chen
# write.table(gd.data_noBatch@dr$pca@cell.embeddings[, 1:2], 'Output/final_seurat_PCA_coords.txt', sep = "\t", quote = F)
# write.table(tpms_without_outliers, 'Output/final_seurat_TPMs_without_outliers_logtransformed.txt', sep = "\t", quote = F)


save.image("./RData/SeuratAnalyses.RData")

