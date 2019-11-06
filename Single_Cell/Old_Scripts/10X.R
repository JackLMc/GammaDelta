library(Seurat)
library(dplyr)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/jlm650/Downloads/filtered_gene_bc_matrices/hg19/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size


sparse.size <- object.size(x = pbmc.data)
sparse.size


dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)


pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)


pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)


pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)


PrintFindClustersParams(object = pbmc)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc)


pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

clusters <- pbmc@ident %>% as.data.frame()
colnames(clusters)[which(names(clusters) == ".")] <- "Cell.Type"
head(clusters)

library(tidyverse)
levels(clusters$Cell.Type)
CD8 <- droplevels(subset(clusters, Cell.Type == "CD8 T cells")) %>% rownames()

# Combining together
pbmc2 <- SubsetData(pbmc, cells.use = CD8, do.clean = T, subset.raw = T)

tpms <- read.table("Data/clean_quant_table.csv", sep = ",", header = T, row.names = 1)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(tpms)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","external_gene_name"),
                values = genes, mart = mart)

row.names(G_list) <- G_list$ensembl_gene_id
G_list <- G_list[, 2, drop = F] 
tpms <- tpms[row.names(G_list), ]
rescale <- function(x) {
  s = sum(x)
  (x/s) * 1e6
}
tpms <- apply(tpms, 2, rescale)
library(plyr)
tpms_with_names <- merge(tpms, G_list, by = 0)
tpms_with_names <- subset(tpms_with_names, select = -c(Row.names))
tpms_with_names <- ddply(tpms_with_names, "external_gene_name", numcolwise(sum))
rownames(tpms_with_names) <- tpms_with_names$external_gene_name
tpms <- subset(tpms_with_names, select = -c(external_gene_name))
#tpms <- log(tpms + 1)


ctrl <- SubsetData(pbmc, cells.use = CD8, do.clean = T, subset.raw = T)
ctrl@meta.data$orig.data <- "PBMC"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)

# Set up stimulated object
stim <- CreateSeuratObject(raw.data = tpms, project = "GAMMADELTA", min.cells = 3, min.genes = 200)
stim@meta.data$orig.data <- "GD"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)


# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))


immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "orig.data", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "orig.data", 
              do.return = TRUE)
plot_grid(p1, p2)


PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


p3 <- MetageneBicorPlot(immune.combined, grouping.var = "orig.data", dims.eval = 1:30, 
                        display.progress = FALSE)



DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)





immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "orig.data", 
                                 dims.align = 1:20)


p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "orig.data", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "orig.data", 
              do.return = TRUE)
plot_grid(p1, p2)






# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = 0.8, dims.use = 1:20)
# Visualization
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "orig.data")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)



