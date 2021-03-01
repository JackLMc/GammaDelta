# A script to investigate the single-cell RNA sequencing data from Gutierrez study
# Author: Jack McMurray
set.seed(2021) # Year of analysis

# Load up packages and make common variables
library(Seurat)
library(tidyverse)

# Read in data
load("./Single_Cell/Output/Projections/Nature_Comms_Stuff.RData")
pbmc.data <- read.delim("./Single_Cell/Data/Nature_Comms/GSE124731_single_cell_rnaseq_gene_counts.txt")
meta <- read.delim("./Single_Cell/Data/Nature_Comms/GSE124731_single_cell_rnaseq_meta_data.txt")
Vd1_meta <- droplevels(subset(meta, cell.type == "Vd1"))

# Vd2_meta <- droplevels(subset(meta, cell.type == "Vd2"))
# Cd8_meta <- droplevels(subset(meta, cell.type == "CD8"))
# Cd4_meta <- droplevels(subset(meta, cell.type == "CD4"))
# iNKT_meta <- droplevels(subset(meta, cell.type == "iNKT"))
# NK_meta <- droplevels(subset(meta, cell.type == "NK"))

pbmc.data <- column_to_rownames(pbmc.data, var = "Gene")
VD1_cells <- pbmc.data[, colnames(pbmc.data) %in% Vd1_meta$cell_id]


library(GenomicFeatures)
# rtracklayer::ucscGenomes()[ , "db"]
# install.packages("RMariaDB")
hg38.ens <- makeTxDbFromUCSC(genome = "hg38", tablename = "knownGene")
exonic <- exonsBy(hg38.ens, by = "gene")
red.exonic <- reduce(exonic)
exon.lengths <- sum(width(red.exonic))
exon.lengths <- exon.lengths %>% as.data.frame() %>% rownames_to_column(., var = "Gene")
colnames(exon.lengths) <- c("Gene", "Length")

# Need to merge in the lengths of the genes...
VD1_cells <- VD1_cells %>% rownames_to_column(., var = "Gene")

head(VD1_cells)


VD1_cells_bind <- merge(VD1_cells, exon.lengths, by = "Gene") %>% column_to_rownames(., var = "Gene")
VD1_tpms <- VD1_cells_bind[!'%in%'(colnames(VD1_cells_bind), c("Length"))] / VD1_cells_bind$Length

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- head(VD1_tpms)
G_list1 <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"),
                values = genes, mart = mart)
VD1_tpms <- VD1_tpms[row.names(G_list1), ]
row.names(G_list1) <- G_list1$ensembl_gene_id
G_list1 <- G_list1[, 2, drop = F] 

rescale <- function(x) {
  s = sum(x)
  (x/s) * 1e6
}

str(VD1_tpms)
VD1_tpms <- apply(VD1_tpms, 2, rescale)

library(plyr)
VD1_with_names <- merge(VD1_tpms, G_list1, by = 0)

VD1_with_names <- subset(VD1_with_names, select = -c(Row.names))
VD1_with_names <- ddply(VD1_with_names, "external_gene_name", numcolwise(sum))
rownames(VD1_with_names) <- VD1_with_names$external_gene_name
VD1_tpms <- subset(VD1_with_names, select = -c(external_gene_name))

tpms <- read.table("./Single_Cell/Data/clean_quant_table.csv", sep = ",", header = T, row.names = 1)

genes <- row.names(tpms)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"),
                 values = genes, mart = mart)
row.names(G_list) <- G_list$ensembl_gene_id
G_list <- G_list[, 2, drop = F] 


tpms <- tpms[row.names(G_list), ]

tpms <- apply(tpms, 2, rescale)
tpms_with_names <- merge(tpms, G_list, by = 0)

tpms_with_names <- subset(tpms_with_names, select = -c(Row.names))
tpms_with_names <- ddply(tpms_with_names, "external_gene_name", numcolwise(sum))
rownames(tpms_with_names) <- tpms_with_names$external_gene_name
tpms <- subset(tpms_with_names, select = -c(external_gene_name))

save.image("./Single_Cell/Output/Projections/Nature_Comms_Stuff.RData")



# Set up control object
NatComm <- CreateSeuratObject(raw.data = VD1_tpms, project = "NatCom", min.cells = 5)
NatComm@meta.data$stim <- "NATURE"
NatComm <- FilterCells(NatComm, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
NatComm <- NormalizeData(NatComm)
NatComm <- ScaleData(NatComm, display.progress = F)
NatComm <- FindVariableGenes(object = NatComm, mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             x.low.cutoff = 2, y.cutoff = 0.1, do.plot = F)

NatComm <- RunPCA(NatComm, do.print = F)
PCAPlot(NatComm, 1, 2)
rotations <- as.data.frame(NatComm@dr$pca@cell.embeddings)
outlier_cells <- rotations[rotations$PC2<(-10), ]
VD1_without_outliers <- VD1_cell[, !colnames(VD1_cell) %in% rownames(outlier_cells)]

NatComm <- CreateSeuratObject(raw.data = VD1_without_outliers,
                              min.cells = 5,project = "NatCom")
NatComm@meta.data$stim <- "NATURE"
NatComm <- NormalizeData(NatComm)
NatComm <- ScaleData(object = NatComm, display.progress = F)
NatComm <- FindVariableGenes(object = NatComm, mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             x.low.cutoff = 2, y.cutoff = 0.1, do.plot = F)

NatComm <- RunPCA(NatComm, do.print = F)
PCAPlot(NatComm, 1, 2)


# Set up stimulated object
GamDelt <- CreateSeuratObject(raw.data = tpms, project = "GAMMADELTA", min.cells = 5)

# GamDelt@meta.data$stim <- "GD"
GamDelt <- FilterCells(GamDelt, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
GamDelt <- NormalizeData(GamDelt)
GamDelt <- ScaleData(GamDelt, display.progress = F)

GamDelt <- FindVariableGenes(object = GamDelt, mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             x.low.cutoff = 2, y.cutoff = 0.1, do.plot = F)

GamDelt <- RunPCA(GamDelt, do.print = F)
PCAPlot(GamDelt, 1, 2)
rotations <- as.data.frame(GamDelt@dr$pca@cell.embeddings)
outlier_cells <- rotations[(rotations$PC1>(20) | rotations$PC2<(-20)), ] # No Outliers
tpms_without_outliers <- tpms[, !colnames(tpms) %in% rownames(outlier_cells)]

# Gene selection for input to CCA
NatComm <- FindVariableGenes(object = NatComm, mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             x.low.cutoff = 2, y.cutoff = 0.1, do.plot = F)

GamDelt <- FindVariableGenes(object = GamDelt, mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             x.low.cutoff = 2, y.cutoff = 0.1, do.plot = F)


g.1 <- head(rownames(NatComm@hvg.info), 1000)
g.2 <- head(rownames(GamDelt@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(NatComm@scale.data))
genes.use <- intersect(genes.use, rownames(GamDelt@scale.data))

name_map <- read.csv("./Single_Cell/Data/name_map.txt", sep = "\t", row.names = 1)
rownames(name_map) = paste('X', rownames(name_map), sep = "")
name_map$sort_class <- ifelse((name_map$sort_class == "naive"), "VD1.CD27HI", "VD1.CD27LO")
head(name_map)

GamDelt <- AddMetaData(object = GamDelt,
                               metadata = name_map,
                               col.name = "stim")

colnames(GamDelt@meta.data)[colnames(GamDelt@meta.data) == "sort_class"] <- "stim"


immune.combined <- RunCCA(NatComm, GamDelt, genes.use = genes.use, num.cc = 30)
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot

p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim", 
              pt.size = 1, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

p1



PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


p3 <- MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:30, 
                        display.progress = FALSE)



DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)


## On full data
immune.combined <- RunPCA(immune.combined, do.print = F)
PCAPlot(immune.combined, 1, 2)
immune.combined <- FindClusters(immune.combined,
                                dims.use = 1:2,
                                resolution = 0.2,
                                print.output = 0,
                                save.SNN = T)

immune.combined@meta.data$stim
?PCAPlot
?DimPlot
PCAPlot(immune.combined, 1, 2, do.return = T, group.by = "stim") +
  theme(legend.position = "top") + ggtitle("Clustered_Class")

DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
                                 dims.align = 1:20)



p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)






# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)
# Visualization
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)




NatComm <- ProjectPCA(NatComm, do.print = F)
PCAPlot(NatComm)
PCHeatmap(NatComm, pc.use = 1, use.full = F,
          do.balanced = T, remove.key = F, col.use = PurpleAndYellow())
NatComm <- FindClusters(NatComm,
                                dims.use = 1:2,
                                resolution = 0.2,
                                print.output = 0,
                                save.SNN = T)

FeaturePlot(object = NatComm, features.plot = c("LTB", "CST7", "GZMB"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), 
            pt.size = 0.5, reduction.use = "pca")


clustered_plot <- PCAPlot(NatComm, 1, 2, do.return = T) 

dat <- as.matrix(NatComm@scale.data)

library(tidyverse)
meta <- as.data.frame(NatComm@meta.data) %>%
  rownames_to_column("Cells")
head(meta)
meta1 <- meta[, c("Cells", "res.0.2")]
meta1$`res.0.2` <- as.factor(meta1$`res.0.2`)
col.cell1 <- c("#009E73","#999999")[meta1$`res.0.2`]
data.frame(meta1$`res.0.2`, col.cell1)



#### Gain the topGenes
library(gplots)
HeatmapGenes <- DimTopGenes(NatComm, dim.use = 1, do.balanced = T, use.full = T)
head(HeatmapGenes)
genes.ordered <- rev(HeatmapGenes)

dat1 <- subset(dat, rownames(dat) %in% HeatmapGenes)
data.use <- as.data.frame(MinMax(data = dat1, min = -2.5, max = 2.5))

#### Better clustering
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

pdf("Figures/Gamma-Delta/Newer/PC1_Heat_plot.pdf", width = 3)
heatmap.2(as.matrix(data.use[genes.ordered,]),
          col = PurpleAndYellow(),
          trace = "none",
          cex.main = 1.5,
          Rowv = F,
          Colv = T,
          ColSideColors = col.cell1,
          scale = "none",
          margin = c(10,5), lhei = c(2,10),
          dendrogram = "none",
          hclustfun = hclustAvg,
          labCol = NA)
par(xpd = T)
legend(x = 0.87, y = 1.065,
       fill = c("#009E73","#999999"),
       legend = levels(meta1$`res.0.2`))
dev.off()










#### Attempting different
#### Use Seurat's clustering method to find clusters in the data.
NatComm <- FindClusters(NatComm,
                        dims.use = 1:2,
                        resolution = 0.2,
                        print.output = 0,
                        save.SNN = T)



NatComm <- RunTSNE(NatComm, reduction.use = "pca", dims.use = 1:5, perplexity=10)
TSNEPlot(NatComm)

