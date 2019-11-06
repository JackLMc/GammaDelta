# T-cell populations
library(Seurat)
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Single_Cell/")

pbmc.data <- Read10X(data.dir = "Data/10X/T_Cell/3k/")
pbmc_3k <- CreateSeuratObject(raw.data = pbmc.data, project = "10X_3K", min.cells = 3, min.genes = 200)
pbmc_3k@meta.data$orig.data <- "3k"

pbmc.data <- Read10X(data.dir = "Data/10X/T_Cell/4k/")
pbmc_4k <- CreateSeuratObject(raw.data = pbmc.data, project = "10X_4K", min.cells = 5)
pbmc_4k@meta.data$orig.data <- "4k"


pbmc.combined <- MergeSeurat(object1 = pbmc_3k, object2 = pbmc_4k, add.cell.id1 = "3K", 
                             add.cell.id2 = "4K", project = "10X_7K")



pmbc.combined <- FilterCells(pbmc.combined, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
pmbc.combined <- NormalizeData(pmbc.combined)


pbmc.combined <- FindVariableGenes(object = pbmc.combined,
                                   mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc.combined@var.genes)
pbmc.combined <- ScaleData(object = pbmc.combined, vars.to.regress = c("nGene", "orig.ident"))

pbmc.combined <- RunPCA(object = pbmc.combined, pc.genes = pbmc.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
                        genes.print = 5)

PrintPCA(object = pbmc.combined, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = pbmc.combined, dim.1 = 1, dim.2 = 2)


pbmc.combined <- FindClusters(object = pbmc.combined, reduction.type = "pca", dims.use = 1:10, 
                              resolution = 1.2, print.output = 0, save.SNN = TRUE, force.recalc = T)

# PrintFindClustersParams(object = pbmc.combined)

pbmc.combined <- RunTSNE(object = pbmc.combined, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc.combined)


pbmc.markers <- FindAllMarkers(object = pbmc.combined, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC) %>% as.data.frame()


cluster1.markers <- FindMarkers(object = pbmc.combined, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(object = pbmc.combined, features.plot = c("CD8A", "CD8B"))
VlnPlot(object = pbmc.combined, features.plot = c("CD4"))


FeaturePlot(object = pbmc.combined, features.plot = c( "CD3E", "PTPRC", "CD8A","CD8B", "CD27", "GZMB", "KLRB1", "CCR7"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
# alphabeta/gammadelta

load("Final_seurat.RData")
str(gd.data_noBatch)
# Gene selection for input to CCA
ctrl <- FindVariableGenes(pbmc.combined, do.plot = F)
stim <- FindVariableGenes(gd.data_noBatch, do.plot = F)
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


