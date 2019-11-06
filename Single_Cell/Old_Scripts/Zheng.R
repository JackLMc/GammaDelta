# T-cell populations
library(Seurat)
setwd("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/4_Gamma_Delta/Single_Cell/")

pbmc.data <- Read10X(data.dir = "Data/Zheng/Cytotoxic/")
pbmc_cy <- CreateSeuratObject(raw.data = pbmc.data, project = "10X_Cyto", min.cells = 5)
pbmc_cy@meta.data$orig.data <- "Cyto"

pbmc.data <- Read10X(data.dir = "Data/Zheng/Naive/")
pbmc_na <- CreateSeuratObject(raw.data = pbmc.data, project = "10X_Naive", min.cells = 5)
pbmc_na@meta.data$orig.data <- "Naive"


pbmc.combined <- MergeSeurat(object1 = pbmc_cy, object2 = pbmc_na, add.cell.id1 = "Cyto", 
                             add.cell.id2 = "Nai", project = "10X_Zheng")



pmbc.combined <- FilterCells(pbmc.combined, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
pmbc.combined <- NormalizeData(pmbc.combined)




pbmc.combined <- FindVariableGenes(object = pbmc.combined,
                                   mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc.combined@var.genes)
pbmc.combined <- ScaleData(object = pbmc.combined, vars.to.regress = c("nUMI"))

pbmc.combined <- RunPCA(object = pbmc.combined, pc.genes = pbmc.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
                        genes.print = 5)

PrintPCA(object = pbmc.combined, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = pbmc.combined, dim.1 = 1, dim.2 = 2)


pbmc.combined <- FindClusters(object = pbmc.combined, reduction.type = "pca", dims.use = 1:10, 
                              resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc = T)

# PrintFindClustersParams(object = pbmc.combined)

pbmc.combined <- RunTSNE(object = pbmc.combined, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc.combined)


pbmc.markers <- FindAllMarkers(object = pbmc.combined, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC) %>% as.data.frame()


cluster1.markers <- FindMarkers(object = pbmc.combined, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(object = pbmc.combined, features.plot = c("CD8A", "CD8B"))
VlnPlot(object = pbmc.combined, features.plot = c("CD4"))


FeaturePlot(object = pbmc.combined, features.plot = c( "CD3E",
                                                       "PTPRC",
                                                       "CD8A",
                                                       "CD8B",
                                                       "CD27",
                                                       "CCR7",
                                                       "GZMB"),
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
# alphabeta/gammadelta

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)