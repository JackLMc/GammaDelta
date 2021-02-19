setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Single_Cell/")
library(Seurat)

## Read in data
tpms <- read.table("Data/clean_quant_table.csv", sep = ",", header = T, row.names = 1)

## Get common names of genes
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(tpms)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","external_gene_name"),
                values = genes, mart = mart)

backup_g_list <- G_list
G_list <- backup_g_list
row.names(G_list) <- G_list$ensembl_gene_id
G_list = G_list[, 2, drop = F] 

### Gain the rows (genes) which have a common name - Genes that are actually genes
dim(tpms)
tpms <- tpms[row.names(G_list), ]
dim(tpms)

## rescaling - Lost gene rows: rescale so they still add to 1 million
rescale <- function(x) {
    s = sum(x)
    (x/s) * 1e6
    }

tpms <- apply(tpms, 2, rescale)

## Merge together
library(plyr)

tpms_with_names <- merge(tpms, G_list, by = 0)
dim(tpms_with_names)

tpms_with_names <- tpms_with_names[, -1]
tpms_with_names <- ddply(tpms_with_names, "external_gene_name", numcolwise(sum))

save.image("Final_seurat.RData")

rownames(tpms_with_names) <- tpms_with_names$external_gene_name

tpms <- tpms_with_names[,-1]

tpms <- log(tpms + 1)

data <- new("seurat", raw.data = tpms)

# check how many genes have at least one transcript in each cell
at_least_one <- apply(tpms, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(tpms),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(tpms, 1, function(x) sum(x>0))
table(tmp>=3)

# all cells have at least 2000 detected genes (Min)
keep <- tmp>=3
tmp <- tpms[keep,]
at_least_one <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one)

# Create object (taking same as Mike)
gd.data <- CreateSeuratObject(raw.data = tpms,
                           min.cells = 3,
                           min.genes = 1000,
                           is.expr = 1,
                           names.field = 2,
                           names.delim = "_",
                           project = "gammadelta")




# Check Mitochondrial genes

# mitochondria genes conveniently start with MT
mito.genes <- grep(pattern = "^MT-", x = rownames(x = gd.data@data), value = TRUE)
length(mito.genes)
percent.mito <- Matrix::colSums(gd.data@raw.data[mito.genes, ]) / Matrix::colSums(gd.data@raw.data)

# check out the meta data
head(gd.data@meta.data)


# add some more meta data
gd.data <- AddMetaData(object = gd.data,
                    metadata = percent.mito,
                    col.name = "percent.mito")

head(gd.data@meta.data)

# plot number of genes, UMIs, and % mitochondria
VlnPlot(object = gd.data,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = gd.data, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = gd.data, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')


table(gd.data@meta.data$percent.mito < 0.05 & gd.data@meta.data$nGene<2500)

#####



data = Setup(data, project = "gammadelta", min.cells = 3, names.field = 2, names.delim = "_", 
             min.genes = 1000, is.expr = 1, do.logNormalize = F) 

data_noBatch = RegressOut(data, latent.vars = c("orig.ident", "nGene"))

data_noBatch=MeanVarPlot(data_noBatch,y.cutoff = 0.1,x.low.cutoff =2 ,fxn.x = expMean,fxn.y = logVarDivMean)

length(data_noBatch@var.genes)

data_noBatch=PCA(data_noBatch,do.print=FALSE)

PCAPlot(data_noBatch, 1, 2)

#get the outliers by inspection
data_noBatch@pca.rot[(data_noBatch@pca.rot$PC2<(-15)|data_noBatch@pca.rot$PC1>20),1:2]

outlier_cells = row.names(data_noBatch@pca.rot[(data_noBatch@pca.rot$PC2<(-15)|data_noBatch@pca.rot$PC1>20),1:2])

tpms_without_outliers = tpms[,!colnames(tpms) %in% outlier_cells]

dim(tpms)

dim(tpms_without_outliers)

data = new("seurat", raw.data=tpms_without_outliers)

data = Setup(data,project="gammadelta",min.cells = 3,names.field = 2,names.delim = "_",min.genes = 1000,is.expr=1,do.logNormalize = F) 

data_noBatch = RegressOut(data, latent.vars=c("orig.ident", "nGene"))

data_noBatch=MeanVarPlot(data_noBatch,y.cutoff = 0.1,x.low.cutoff =2 ,fxn.x = expMean,fxn.y = logVarDivMean)

length(data_noBatch@var.genes)

data_noBatch=PCA(data_noBatch,do.print=FALSE)

PCAPlot(data_noBatch, 1, 2)

VizPCA(data_noBatch,1:2)

PCHeatmap(data_noBatch,pc.use = 1:2,do.balanced = TRUE)

data_noBatch=JackStraw(data_noBatch,num.replicate = 100,do.print = FALSE) 

JackStrawPlot(data_noBatch,PCs = 1:10)

PCElbowPlot(data_noBatch)

data_noBatch=ProjectPCA(data_noBatch,do.print=FALSE)

PCHeatmap(data_noBatch,pc.use = 1:2,use.full = TRUE,do.balanced = TRUE,remove.key = TRUE)

name_map = read.csv("../name_map.txt", sep="\t", row.names=1)

head(name_map)

rownames(name_map) = paste('X', rownames(name_map), sep="")

name_map

head(name_map)

data_noBatch <- AddMetaData(data_noBatch, name_map[rownames(data_noBatch@data.info),'sort_class',drop=FALSE], "sort_class")

head(data_noBatch@data.info)

PCAPlot(data_noBatch, 1, 2, group.by="sort_class")

FeaturePlot(data_noBatch, c('PRF1', 'GZMA', 'CX3CR1', 'SELL', 'TCF7', 'LEF1', 'IL7R', 'LTB'), reduction.use='pca')

data_noBatch = FindClusters(data_noBatch, pc.use=1:2, resolution=0.2, print.output = 0, save.SNN = T)
PCAPlot(data_noBatch, 1, 2)

current.cluster.ids = c(0,1)
cluster.ids = c("EMRA_clust", "naive_clust")

data_noBatch@ident <- plyr::mapvalues(data_noBatch@ident, from = current.cluster.ids, to = cluster.ids)

plot1 <- PCAPlot(data_noBatch, do.return=T, group.by="sort_class", no.legend=T, cols.use = c("#e41a1c", "#377eb8"), do.label=T)
plot2 <- PCAPlot(data_noBatch, do.return=T, no.legend=T, cols.use = c("#e41a1c", "#377eb8"), do.label=T)
MultiPlotList(list(plot1, plot2), cols=2)

weightings = data_noBatch@pca.x

PC1_orders = weightings[order(weightings$PC1),'PC1',drop=F]

PC2_orders = weightings[order(weightings$PC2),'PC2',drop=F]

write.table(PC1_orders, 'final_seurat_PC1_weightings.txt', sep="\t", quote=F)

write.table(PC2_orders, 'final_seurat_PC2_weightings.txt', sep="\t", quote=F)

EMRA_cluster.markers <- FindMarkers(data_noBatch, ident.1 = 'EMRA_clust', thresh.use = 0.25, only.pos = T)

naive_cluster.markers <- FindMarkers(data_noBatch, ident.1 = 'naive_clust', thresh.use = 0.25, only.pos = T)

write.table(EMRA_cluster.markers,'final_seurat_EMRA_cluster_markers.txt', sep="\t", quote=F)
write.table(naive_cluster.markers,'final_seurat_naive_cluster_markers.txt', sep="\t", quote=F)

?write.table

write.table(data_noBatch@ident,'final_seurat_cluster_assignments.txt', sep="\t", quote=F)

write.table(data_noBatch@pca.rot[,1:2], 'final_seurat_PCA_coords.txt', sep="\t", quote=F)

write.table(tpms_without_outliers, 'final_seurat_TPMs_without_outliers_logtransformed.txt', sep="\t", quote=F)

save.image("Final_seurat.RData")

load("Final_seurat.RData")

slotNames(data_noBatch)

head(data_noBatch@data)

FeaturePlot(data_noBatch, c('MKI67'), reduction.use = 'pca')

hist(as.numeric(tpms['MKI67',]))
