# A script to investigate the single-cell RNA sequencing data from Mike Stubbington on GammaDelta T-cells
# Original Author: Xi Chen (pre Seurat Version 2.3.4)#
# Co-Author: Jack McMurray (Updated to Seurat 2.3.4, and then to version 4.0.0#

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
# load("RData/SeuratAnalyses.RData")

#### Read in the table of TPM expression values for the cells that passed QC - 
##### This was performed outside of R
tpms <- read.table("Data/clean_quant_table.csv", sep = ",", header = T, row.names = 1)

#### Rename the gene names from ENSEMBL identifiers to gene symbols to make it easier to interpret results
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(tpms)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"),
                values = genes, mart = mart)

backup_g_list <- G_list
G_list <- backup_g_list
row.names(G_list) <- G_list$ensembl_gene_id
G_list <- G_list[, 2, drop = F] 


dim(tpms)
tpms <- tpms[row.names(G_list), ]
dim(tpms)

#### We now have a TPM table with gene names instead of ENSEMBL IDs. 
#### We've lost 1434 rows because they're strange entries which aren't really genes.
#### This is nearly double what was lost in original analysis (620)
#### We can safely ignore these but it's better to rescale the columns so they sum to 1 million.
rescale <- function(x) {
    s = sum(x)
    (x/s) * 1e6}

tpms <- apply(tpms, 2, rescale)

#### There are cases where more than one ENSEMBL ID resolves to the same gene name
#### We want to add the TPMs for each duplicated name.
library(plyr)

tpms_with_names <- merge(tpms, G_list, by = 0)
tpms_with_names <- subset(tpms_with_names, select = -c(Row.names))
tpms_with_names <- plyr:: ddply(tpms_with_names, "external_gene_name", numcolwise(sum))

rownames(tpms_with_names) <- tpms_with_names$external_gene_name
tpms <- subset(tpms_with_names, select = -c(external_gene_name))

## Setting up Seurat object ##
#### Normalisation is for Counts, so avoid doing that and pass in log transformed TPMs
tpms <- log(tpms + 1)

#### Check how many genes have at least one transcript in each cell
at_least_one <- apply(tpms, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(tpms),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

#### Manually check the number of genes detected in three or more cells
#### A lot of genes are not detected in 3 or more cells
tmp <- apply(tpms, 1, function(x) sum(x>0))
table(tmp>=3)
#### all cells have at least 2000 detected genes (Min)

keep <- tmp>=3
tmp <- tpms[keep, ]
at_least_one <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one)

#### Create the object with new Seurat workflow
#### Using the parameters that Mike used in his script (not informed from above)
#### Can't be done using `Setup` anymore, function is deprecated
gd.data <- CreateSeuratObject(counts = tpms,
                              min.cells = 3,
                              min.genes = 1000,
                              is.expr = 1,
                              names.field = 2,
                              names.delim = "_",
                              project = "gammadelta")

## Check Mitochondrial genes ##
#### A good QC metric is the percentage of mitochrondrial RNA
#### Cells that were lysed, will lose cytoplasmic RNA, but mito will be kept and sequenced
gd.data[["percent.mt"]] <- PercentageFeatureSet(gd.data, pattern = "^MT-")

VlnPlot(gd.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gd.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

gd.data <- subset(gd.data, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5) # Choose these based on previous plots

#### Not a problem, all cells have below 5% mitochondrial RNA

gd.data <- FindVariableFeatures(gd.data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gd.data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gd.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
plot1 + plot2


#### Although, initially, I thought that there weren't batch effects in the data.
#### It turns out that there are and that Donor 6 sits apart from the others.
#### This makes sense since that plate was done at a very different time.
#### So, I want to remove the batch effects (this is nice and easy because the plates aren't confounded).
#### Exploratory analyses found that the number of genes detected influences the results and that the biology appears easier to interpret if we remove this.
all.genes <- rownames(gd.data)
gd.data <- ScaleData(gd.data, features = all.genes, vars.to.regress = c("orig.ident", "nCount_RNA"))
# Regressing out the orig.ident and the nCount_RNA as to remove any batch effects

# SKIP THE NORMALISATION STEP HERE. SINCE I PASSED IT LOG(TPM + 1)

#### Next we want to find a set of genes that are sufficiently variable
#### (and expressed at sufficient levels) to be interesting.
#### This is a *bit* arbritary but these cutoffs are fairly relaxed so I don't think we're missing much.
#### Making them more relaxed tended to make the results look weird and not in an interesting way.

## Initial analyses ##
#### Run PCA
gd.data <- RunPCA(gd.data, features = VariableFeatures(object = gd.data))
print(gd.data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(gd.data, dims = 1:2, reduction = "pca")
DimPlot(gd.data, reduction = "pca")

### Produce some heatmaps about the data
DimHeatmap(gd.data, dims = 1, cells = 445, balanced = T)

#### Find Dimensionality of the data
gd.data <- JackStraw(gd.data, num.replicate = 100)
gd.data <- ScoreJackStraw(gd.data, dims = 1:20)

JackStrawPlot(gd.data, dims = 1:15)

# PC1 and PC12
DimPlot(gd.data, reduction = "pca", dims = c(1, 12))


gd.data <- FindNeighbors(gd.data, dims = 1:10)
gd.data <- FindClusters(gd.data, resolution = 0.5)


## Run umad
# reticulate::py_install(packages = 'umap-learn')
gd.data <- RunUMAP(gd.data, dims = 1:10)
DimPlot(gd.data, reduction = "umap") # Doen't really improve it

save.image("RData/SeuratAnalyses.RData")





# #### How many genes does this give us?
# length(gd.data_noBatch@var.genes)
# #### Less than what was said on the original analysis (4,573)
# 
# #### Perform PCA using these genes to see what the populations look like
# gd.data_noBatch <- RunPCA(gd.data_noBatch, do.print = F)
# 
# # Find Outlier cells
# PCAPlot(gd.data_noBatch, 1, 2)
# #### There are outliers, which make it difficult to see the strucutre in the main population.
# #### These are likely technical artefacts, so we'll remove them.
# 
# rotations <- as.data.frame(gd.data_noBatch@dr$pca@cell.embeddings)
# outlier_cells <- rotations[(rotations$PC1>(20) | rotations$PC2<(-20)), ]
# tpms_without_outliers <- tpms[, !colnames(tpms) %in% rownames(outlier_cells)]
# 
# dim(tpms)
# dim(tpms_without_outliers)
# ## The 3 outlier cells have been removed
# 
# ## Repeating ##
# #### Doing everything above, but without the outlier cells
# gd.data <- CreateSeuratObject(raw.data = tpms_without_outliers,
#                               min.cells = 3,
#                               min.genes = 1000,
#                               is.expr = 1,
#                               names.field = 2,
#                               names.delim = "_",
#                               project = "gammadelta")
# 
# gd.data_noBatch <- ScaleData(object = gd.data,
#                              vars.to.regress = c("orig.ident", "nGene"))
# 
# gd.data_noBatch <- FindVariableGenes(object = gd.data_noBatch,
#                                      mean.function = ExpMean,
#                                      dispersion.function = LogVMR,
#                                      x.low.cutoff = 2,
#                                      y.cutoff = 0.1)
# 
# 
# length(gd.data_noBatch@var.genes)
# 
# #### Perform the PCA again, without outliers
# gd.data_noBatch <- RunPCA(gd.data_noBatch, do.print = F)
# PCAPlot(gd.data_noBatch, 1, 2)
# 
# #### We can have a look at the genes with the heighest weights in the first two PC's to see what bio they represent
# VizPCA(gd.data_noBatch, 1:2)
# 
# PCHeatmap(gd.data_noBatch, 
#           pc.use = 1:2, 
#           do.balanced = T)
# 
# #### It looks like PC1 describes most bio that we're expecting to see.
# #### PC2 is less clear (we'll come back to this)
# 
# #### Next we can try and determine how many PCs contain use information about our data
# gd.data_noBatch <- JackStraw(gd.data_noBatch,
#                              num.replicate = 100) 
# pdf("Figures/Gamma-Delta/JackStrawPlot.pdf")
# JackStrawPlot(gd.data_noBatch,
#               PCs = 1:6)
# dev.off()
# 
# pdf("Figures/Gamma-Delta/ElbowPlot.pdf")
# PCElbowPlot(gd.data_noBatch)
# dev.off()
# 
# #### So, it looks like only PC1 is of interest. This means that we can ignore more complicated dimension-reduction techniques
# #### Like tSNE and focus on PC1 (maybe PC2)
# 
# #### We calculated the PCAs using a subset of seemingly variable genes.
# #### We now project the PCA analysis onto the entire dataset to give us weightings for all the genes.
# #### This helps to ensure that we're not missing anything by filtering the genes
# gd.data_noBatch <- ProjectPCA(gd.data_noBatch, do.print = F)
# 
# getwd()
# pdf("Figures/Gamma-Delta/PC_Contrib.pdf")
# PCHeatmap(gd.data_noBatch, pc.use = 1, use.full = T,
#           do.balanced = T, remove.key = F, col.use = PurpleAndYellow())
# dev.off()
# 
# 
# #### The PC1 heatmap is really nice. This is pretty much as we saw before
# #### PC2 isn't particularly interesting
# 
# ## Interpretin the PCA ##
# #### We want to check that the separation of cells on PC1 makes sense.
# #### Lets see how the sort idenities (CD27+ or CD27-) map onto the plot
# name_map <- read.csv("Data/name_map.txt", sep = "\t", row.names = 1)
# rownames(name_map) = paste('X', rownames(name_map), sep = "")
# name_map$sort_class <- ifelse((name_map$sort_class == "naive"), "VD1.CD27HI", "VD1.CD27LO")
# head(name_map)
# 
# gd.data_noBatch <- AddMetaData(object = gd.data_noBatch,
#             metadata = name_map,
#             col.name = "sort_class")
# 
# 
# PC <- PCAPlot(gd.data_noBatch, 1, 2, group.by = "sort_class", do.return = T, cols.use = cbcols) +
#   theme(legend.position = "top") + ggtitle("Sorted_Class")
# ggsave("Figures/Gamma-Delta/Sorted_PCA.pdf",
#        plot = PC,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# #### Looks pretty convincing that PC1 explains the difference between the two original populations.
# #### Use Seurat's clustering method to find clusters in the data.
# gd.data_noBatch <- FindClusters(gd.data_noBatch,
#                              dims.use = 1:2,
#                              resolution = 0.2,
#                              print.output = 0,
#                              save.SNN = T)
# 
# #### Do these split along the same lines as the sort classes
# current.cluster.ids <- c(0, 1)
# cluster.ids <- c("Cluster_1", "Cluster_2")
# 
# gd.data_noBatch@ident <- plyr::mapvalues(gd.data_noBatch@ident,
#                                          from = current.cluster.ids,
#                                          to = cluster.ids)
# 
# clustered_plot <- PCAPlot(gd.data_noBatch, 1, 2, do.return = T, cols.use = cbcols) +
#   theme(legend.position = "top") + ggtitle("Clustered_Class")
# 
# #### NEWER ONES FOR EACH!
# embeddings.use <- GetDimReduction(object = gd.data_noBatch, reduction.type = "pca", slot = "cell.embeddings")
# data.plot <- as.data.frame(x = embeddings.use)
# 
# ### Clustered
# clusters <- as.data.frame(gd.data_noBatch@ident)
# data.plot_cluster <- merge(data.plot, clusters, by = "row.names", all.x = T)
# colnames(data.plot_cluster)[which(names(data.plot_cluster) == "gd.data_noBatch@ident")] <- "Clustered_Class"
# 
# 
# cols_for_clusts <- c("Cluster_2" = "#0072B2",
#                      "Cluster_1" = "#D55E00")
# 
# triangle_sort <- ggplot(data.plot_cluster, mapping = aes(x = PC1, y = PC2, colour = Clustered_Class)) + 
#   geom_point(shape = 17) + 
#   scale_color_manual(values = cols_for_clusts) + 
#   theme(legend.position = "top") +
#   theme(legend.title = element_blank())
# 
# ggsave("Figures/Gamma-Delta/Newer/Triangle_Clustered.pdf",
#        plot = triangle_sort,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# ### Sorted
# clusters <- as.data.frame(gd.data_noBatch@meta.data)
# data.plot_sort <- merge(data.plot, clusters, by = "row.names", all.x = T)
# 
# circl_sort <- ggplot(data.plot_sort, mapping = aes(x = PC1, y = PC2, colour = sort_class)) + 
#   geom_point(shape = 16) + 
#   scale_color_manual(values = cbcols) + 
#   theme(legend.position = "top") +
#   theme(legend.title = element_blank())
# 
# ggsave("Figures/Gamma-Delta/Newer/Circle_Sorted.pdf",
#        plot = circl_sort,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# 
# 
# #### Plot old ones together
# plot_grid(PCAPlot(gd.data_noBatch,
#                   do.return = T,
#                   group.by = "sort_class",
#                   no.legend = T,
#                   cols.use = cbcols,
#                   do.label = T),
#           PCAPlot(gd.data_noBatch,
#                   do.return = T,
#                   no.legend = T,
#                   cols.use = cbcols,
#                   do.label = T))
# 
# #### Find the overlap between them
# # library(plotly)
# # p <- PCAPlot(gd.data_noBatch, 1, 2, do.return = T, do.hover = T, data.hover = "sort_class")
# # 
# # setwd("./Figures")
# # htmlwidgets::saveWidget(as_widget(p), "PCA_hover.html")
# # setwd("../")
# # 
# # Sys.setenv("plotly_username" = "JackMcMurray")
# # Sys.setenv("plotly_api_key" = "nFWGB59b0FzSipM16Sgt")
# # api_create(x = p, filename = "Gamma_Delta")
# 
# Overlap <- PCAPlot(gd.data_noBatch, 1, 2, pt.shape = "sort_class", do.return = T) +
#   theme(legend.position = "top")
# ggsave("Figures/Gamma-Delta/Overlapping_PCA.png",
#        plot = Overlap,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# #### *Yes!*
# #### There are some in the centre that aren't concordant, but this is either:
# #### A limitation of the clustering
# #### Representing the cells that were borderline of the FACS gates
# #### The index sorting data would be really useful here to tease this apart
# library(tidyverse)
# clustered <- gd.data_noBatch@ident %>% as.data.frame() %>% rownames_to_column(var = "cell.ID")
# colnames(clustered)[which(names(clustered) == ".")] <- "clustered_class"
# clustered$matches_ <- ifelse((clustered$clustered_class == "Cluster_2"), "VD1.CD27LO", "VD1.CD27HI")
# 
# sorted <- gd.data_noBatch@meta.data %>% as.data.frame() %>% rownames_to_column(var = "cell.ID")
# sorted1 <- sorted[, c("cell.ID", "sort_class")]
# 
# comparison <- merge(clustered, sorted1, by = "cell.ID")
# comparison$clustered_class <- as.character(comparison$clustered_class)
# 
# total.cells <- length(comparison$cell.ID)
# 
# matching <- subset(comparison, sort_class == matches_)
# matching.cells <- length(matching$cell.ID)
# 
# concordance <- matching.cells/total.cells
# concordance
# 
# ## Remake the Heatmap ##
# #### Get scaled data, and the data to colour by
# dat <- as.matrix(gd.data_noBatch@scale.data)
# 
# library(tidyverse)
# meta <- as.data.frame(gd.data_noBatch@meta.data) %>%
#   rownames_to_column("Cells")
# head(meta)
# meta1 <- meta[, c("Cells", "sort_class")]
# meta1$sort_class <- as.factor(meta1$sort_class)
# col.cell1 <- c("#009E73","#999999")[meta1$sort_class]
# data.frame(meta1$sort_class, col.cell1)
# 
# 
# 
# #### Gain the topGenes
# library(gplots)
# HeatmapGenes <- DimTopGenes(gd.data_noBatch, dim.use = 1, do.balanced = T, use.full = T)
# head(HeatmapGenes)
# genes.ordered <- rev(HeatmapGenes)
# 
# dat1 <- subset(dat, rownames(dat) %in% HeatmapGenes)
# data.use <- as.data.frame(MinMax(data = dat1, min = -2.5, max = 2.5))
# 
# #### Better clustering
# distCor <- function(x) as.dist(1-cor(t(x)))
# hclustAvg <- function(x) hclust(x, method = "average")
# 
# pdf("Figures/Gamma-Delta/Newer/PC1_Heat_plot.pdf", width = 3)
# heatmap.2(as.matrix(data.use[genes.ordered,]),
#           col = PurpleAndYellow(),
#           trace = "none",
#           cex.main = 1.5,
#           Rowv = F,
#           Colv = T,
#           ColSideColors = col.cell1,
#           scale = "none",
#           margin = c(10,5), lhei = c(2,10),
#           dendrogram = "none",
#           hclustfun = hclustAvg,
#           labCol = NA)
# par(xpd = T)
# legend(x = 0.87, y = 1.065,
#        fill = c("#009E73","#999999"),
#        legend = levels(meta1$sort_class))
# dev.off()
# 
# ## Marker genes ##
# #### Since PC1 appears to contain all the biology of interest, the best way to look for genes that are of interest
# #### is to look at the PC1 weightings for the genes
# #### Very negative weightings = EMRA-like cells, positive weightings = naive cells
# #### How do some markers of interest look
# pdf("Figures/Gamma-Delta/Newer/Mark_of_Interest.pdf", width = 18)
# FeaturePlot(gd.data_noBatch, c("NKG7",
#                                "GZMB",
#                                "PRF1",
#                                "GNLY",
#                                "GZMA",
#                                "CX3CR1",
#                                "LTB",
#                                "TCF7",
#                                "LEF1",
#                                "IL7R",
#                                "NOSIP",
#                                "SELL"), 
#             reduction.use = "pca", nCol = 6)
# dev.off()
# 
# #### As expected from the heatmaps, these look good
# 
# #### I looked at the genes with strong weightings in either direction in PC2 and they looked metabolic
# weightings <- as.data.frame(gd.data_noBatch@dr$pca@gene.loadings.full)
# 
# PC1_orders <- weightings[order(weightings$PC1), "PC1", drop = F]
# PC2_orders <- weightings[order(weightings$PC2), "PC2", drop = F]
# 
# write.table(PC1_orders, "Output/final_seurat_PC1_weightings.txt", sep = "\t", quote = F)
# write.table(PC2_orders, "Output/final_seurat_PC2_weightings.txt", sep = "\t", quote = F)
# 
# sdev <- gd.data_noBatch@dr$pca@sdev
# 
# ### Find the correlataions
# #### Functions
# var_cor_func <- function(var.loadings, comp.sdev){
#   var.loadings*comp.sdev
# }
# contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
# 
# var.coord <- t(apply(weightings, 1, var_cor_func, sdev))
# head(var.coord)
# ## Calculate the Cos2 (square of the coordinates)
# var.cos2 <- var.coord^2
# 
# ## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
# comp.cos2 <- apply(var.cos2, 2, sum)
# var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
# contrib.var <- as.data.frame(var.contrib) %>% rownames_to_column(var = "GeneID")
# nrow(contrib.var)
# 
# 
# try <- contrib.var[contrib.var$GeneID %in% HeatmapGenes,]
# sum(try$PC1)
# 
# 
# ## Choosing the right model is important here
# ### I think I need a negative binomial, but it throws errors because I'm not using counts...
# gd.data_noBatch@data <- as.matrix(gd.data_noBatch@data)
# cluster.markers <- FindMarkers(gd.data_noBatch,
#                                ident.1 = WhichCells(gd.data_noBatch, ident = "Cluster_2"),
#                                ident.2 = WhichCells(gd.data_noBatch, ident = "Cluster_1"),
#                                test.use = "bimod",
#                                only.pos = F)
# 
# 
# cluster.markers$Direction <- ifelse((cluster.markers$avg_logFC > 0), "effector/Cluster2", "naive/Cluster1")
# cluster.markers <- cluster.markers %>% rownames_to_column(., var  = "SYMBOL")
# 
# write.csv(cluster.markers, file = "./Output/SC_DE.csv")
# 
# sigs <- droplevels(subset(cluster.markers, p_val_adj <= 0.051))
# 
# 
# # library(devtools)
# # install_github("https://github.com/RGLab/MAST")
# 
# myData <- as.data.frame(cluster.markers) %>% rownames_to_column(., var = "SYMBOL")
# myData$padjThresh <- as.factor(myData$p_val_adj < 0.05)
# myData$Labels <- myData$SYMBOL
# labelled_genes <- c("CX3CR1", "GZMA", "PRF1", "SELL", "IL7R", "CCR7")
# 
# myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
# 
# pdf(file = "./Single_Cell/Figures/Gamma-Delta/Paper_ready/1C_Volcano.pdf", width = 6, height = 6)
# library(ggrepel)
# ggplot(data = myData, aes(x = avg_logFC, y = -log10(p_val))) + 
#   geom_point(alpha = 0.2, size = 1) +
#   theme(legend.position = "none", text = element_text(size = 10)) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ## Colours for significance
#   geom_point(data = subset(myData, padjThresh==T & avg_logFC<(-0.5)), aes(avg_logFC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "blue4") +
#   geom_point(data = subset(myData, padjThresh==T & avg_logFC>0.5), aes(avg_logFC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "orangered3") +
#   ## Label these
#   geom_label_repel(mapping = aes(label = Labels), box.padding = 0.35, point.padding = 0.2) +
#   # geom_text_repel(data=subset(myData,padjThresh==TRUE & logFC>0.5), aes(logFC,-log10(P.Value),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
#   labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
# dev.off()
# 
# myData <- as.data.frame(CD27LO.vs.CD27HI)
# myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)
# myData$Labels <- myData$SYMBOL
# 
# labelled_genes <- c("LTB", "LEF1", "TCF7", "CCR7", "SELL", "CD27",
#                     "GZMB", "PRF1", "TBX21", "PRDM1", "NKG7", "GNLY")
# myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
# 
# 
# 
# 
# 
# 
# 
# #### We can also ask Seurat to define marker genes for the two observed clusters and save those to a file
# write.table(EMRA_cluster.markers, "Output/final_seurat_EMRA_cluster_markers.txt", sep = "\t", quote = F)
# write.table(naive_cluster.markers, "Output/final_seurat_naive_cluster_markers.txt", sep = "\t", quote = F)
# 
# #### Also write out the cluster assignments
# write.table(gd.data_noBatch@ident, 'Output/final_seurat_cluster_assignments.txt', sep = "\t", quote = F)
# 
# 
# #### Write out some tables that I want in other analyses in Python
# write.table(gd.data_noBatch@dr$pca@cell.embeddings[, 1:2], 'Output/final_seurat_PCA_coords.txt', sep = "\t", quote = F)
# write.table(tpms_without_outliers, 'Output/final_seurat_TPMs_without_outliers_logtransformed.txt', sep = "\t", quote = F)
# 
# 
# # save.image("Final_seurat.RData")
# load("./Single_Cell/Final_seurat.RData")
# 
# #### All Plots, One Place ####
# # Figure 1A - Elbow Plot
# pdf("Figures/Gamma-Delta/Paper_ready/1A_ElbowPlot.pdf")
# PCElbowPlot(gd.data_noBatch)
# dev.off()
# 
# # Figure 1B - Clustered class
# embeddings.use <- GetDimReduction(object = gd.data_noBatch, reduction.type = "pca", slot = "cell.embeddings")
# data.plot <- as.data.frame(x = embeddings.use)
# 
# clusters <- as.data.frame(gd.data_noBatch@ident)
# data.plot_cluster <- merge(data.plot, clusters, by = "row.names", all.x = T)
# colnames(data.plot_cluster)[which(names(data.plot_cluster) == "gd.data_noBatch@ident")] <- "Clustered_Class"
# 
# cols_for_clusts <- c("Cluster_2" = "#999999",
#                      "Cluster_1" = "#009E73")
# 
# Circle_Clust <- ggplot(data.plot_cluster, mapping = aes(x = PC1, y = PC2, colour = Clustered_Class)) + 
#   geom_point(shape = 16) + 
#   scale_color_manual(values = cols_for_clusts) + 
#   theme(legend.position = "top") +
#   theme(legend.title = element_blank())
# dev.off()
# ggsave("Figures/Gamma-Delta/Paper_ready/1B_ClusteredCells.pdf",
#        plot = Circle_Clust,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# # 1C - Volcano Plot
# ## Choosing the right model is important here
# ### I think I need a negative binomial, but it throws errors because I'm not using counts...
# gd.data_noBatch@data <- as.matrix(gd.data_noBatch@data)
# cluster.markers <- FindMarkers(gd.data_noBatch,
#                                ident.1 = WhichCells(gd.data_noBatch, ident = "Cluster_2"),
#                                ident.2 = WhichCells(gd.data_noBatch, ident = "Cluster_1"),
#                                resolution = 0.2, print.output = 0,
#                                save.SNN = T,
#                                only.pos = F)
# 
# # library(devtools)
# # install_github("https://github.com/RGLab/MAST")
# myData <- as.data.frame(cluster.markers) %>% 
#   rownames_to_column(., var = "SYMBOL")
# myData$padjThresh <- as.factor(myData$p_val_adj < 0.05)
# 
# Single_cell <- droplevels(subset(myData, p_val_adj < 0.05))[, c("SYMBOL", "avg_logFC", "p_val_adj")]
# Single_cell$Direction <- ifelse((Single_cell$avg_logFC > 0), "effector/Cluster2", "naive/Cluster1")
# write.csv(x = Single_cell, file = "./Single_Cell/Output/SC_DE.csv", row.names = F)
# 
# myData$Labels <- myData$SYMBOL
# labelled_genes <- c("CX3CR1", "GZMA", "PRF1", "SELL", "IL7R", "CCR7", "TCF7", "LEF1", "TBX21", "EOMES")
# 
# myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
# 
# pdf(file = "./Single_Cell/Figures/Gamma-Delta/Paper_ready/1E_Volcano.pdf", width = 8, height = 8)
# library(ggrepel)
# ggplot(data = myData, aes(x = avg_logFC, y = -log10(p_val))) + 
#   geom_point(alpha = 0.2, size = 1) +
#   theme(legend.position = "none", text = element_text(size = 10)) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ## Colours for significance
#   geom_point(data = subset(myData, padjThresh==T & avg_logFC<(-0.5)), aes(avg_logFC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "blue4") +
#   geom_point(data = subset(myData, padjThresh==T & avg_logFC>0.5), aes(avg_logFC, -log10(p_val)), alpha = 0.6, size = 0.6, colour = "orangered3") +
#   ## Label these
#   geom_label_repel(mapping = aes(label = Labels), box.padding = 0.35, point.padding = 0.2) +
#   # geom_text_repel(data=subset(myData,padjThresh==TRUE & logFC>0.5), aes(logFC,-log10(P.Value),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
#   labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
# dev.off()
# 
# 
# ## Think this is from the bulk stuff for help
# # myData <- as.data.frame(CD27LO.vs.CD27HI)
# # myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)
# # myData$Labels <- myData$SYMBOL
# # 
# # labelled_genes <- c("LTB", "LEF1", "TCF7", "CCR7", "SELL", "CD27",
# #                     "GZMB", "PRF1", "TBX21", "PRDM1", "NKG7", "GNLY")
# # myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
# 
# 
# 
# # 1D - Heatmap of top
# dat <- as.matrix(gd.data_noBatch@scale.data)
# 
# library(tidyverse)
# 
# meta <- as.data.frame(gd.data_noBatch@ident) %>%
#   rownames_to_column("Cells")
# head(meta)
# colnames(meta) <- c("Cells", "Clustered_Class")
# meta1 <- meta[, c("Cells", "Clustered_Class")]
# meta1$Clustered_Class <- as.factor(meta1$Clustered_Class)
# col.cell1 <- c("#009E73","#999999")[meta1$Clustered_Class]
# data.frame(meta1$Clustered_Class, col.cell1)
# 
# library(gplots)
# head(HeatmapGenes)
# HeatmapGenes <- DimTopGenes(gd.data_noBatch, dim.use = 1, do.balanced = T, use.full = T)
# genes.ordered <- rev(HeatmapGenes)
# 
# dat1 <- subset(dat, rownames(dat) %in% HeatmapGenes)
# data.use <- as.data.frame(MinMax(data = dat1, min = -2.5, max = 2.5))
# 
# #### Better clustering
# distCor <- function(x) as.dist(1-cor(t(x)))
# hclustAvg <- function(x) hclust(x, method = "average")
# 
# pdf("Figures/Gamma-Delta/Paper_ready/1D_PC1_Heat_plot.pdf", width = 3)
# heatmap.2(as.matrix(data.use[genes.ordered,]),
#           col = PurpleAndYellow(),
#           trace = "none",
#           cex.main = 1.5,
#           Rowv = F,
#           Colv = T,
#           ColSideColors = col.cell1,
#           scale = "none",
#           margin = c(10,5), lhei = c(2,10),
#           dendrogram = "none",
#           hclustfun = hclustAvg,
#           labCol = NA)
# par(xpd = T)
# legend(x = 0.87, y = 1.065,
#        fill = c("#009E73","#999999"),
#        legend = levels(meta1$Clustered_Class))
# dev.off()
# 
# pdf("./Single_Cell/Figures/Gamma-Delta/Paper_ready/1E_Mark_of_Interest.pdf", width = 21)
# FeaturePlot(gd.data_noBatch, c("NKG7", "GZMB", "PRF1", "GNLY", "GZMA", "CX3CR1",
#                                "LTB", "TCF7", "LEF1", "IL7R", "NOSIP", "SELL"), 
#             reduction.use = "pca", nCol = 6)
# dev.off()
# 
# # 1F - Sorted class confirmation
# embeddings.use <- GetDimReduction(object = gd.data_noBatch, reduction.type = "pca", slot = "cell.embeddings")
# data.plot <- as.data.frame(x = embeddings.use)
# 
# clusters <- as.data.frame(gd.data_noBatch@meta.data)
# data.plot_sort <- merge(data.plot, clusters, by = "row.names", all.x = T)
# 
# Triangle_Sort <- ggplot(data.plot_sort, mapping = aes(x = PC1, y = PC2, colour = sort_class)) + 
#   geom_point(shape = 17) + 
#   scale_color_manual(values = cbcols) + 
#   theme(legend.position = "top") +
#   theme(legend.title = element_blank())
# 
# ggsave("Figures/Gamma-Delta/Paper_ready/1F_SortedCells.pdf",
#        plot = Triangle_Sort,
#        dpi = 600,
#        width = 6,
#        height = 6)
# 
# # 1G - TRaCeR
