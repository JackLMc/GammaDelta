setwd("/Users/jlm650/Downloads/filtered_gene_bc_matrices/GRCh38/")
library(Seurat)

data_dir <- '/Users/jlm650/Downloads/filtered_gene_bc_matrices/GRCh38/'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
expression_matrix <- log(expression_matrix + 1)


#### Check how many genes have at least one transcript in each cell
at_least_one <- apply(expression_matrix, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(expression_matrix),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

seurat_object <- CreateSeuratObject(raw.data = expression_matrix,
                                    names.field = 2,
                                    names.delim = "_",
                                    project = "PBMC")


## Check Mitochondrial genes ##
#### A good QC metric is the percentage of mitochrondrial RNA
#### Cells that were lysed, will lose cytoplasmic RNA, but mito will be kept and sequenced
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_object@data), value = T)
length(mito.genes)
percent.mito <- Matrix::colSums(seurat_object@raw.data[mito.genes, ]) / Matrix::colSums(seurat_object@raw.data)

#### Add the percentage mitochondrial DNA to the meta data
seurat_object <- AddMetaData(object = seurat_object,
                       metadata = percent.mito,
                       col.name = "percent.mito")

VlnPlot(object = seurat_object,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = seurat_object, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = seurat_object, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

table(seurat_object@meta.data$percent.mito < 0.05)
seurat_object <- FilterCells(object = seurat_object,
                       subset.names = c("percent.mito"),
                       low.thresholds = c(-Inf),
                       high.thresholds = c(0.05))

str(seurat_object)
seurat_object_noBatch <- ScaleData(object = seurat_object,
                             vars.to.regress = c("orig.ident", "nGene"))


#### Next we want to find a set of genes that are sufficiently variable
#### (and expressed at sufficient levels) to be interesting.
#### This is a *bit* arbritary but these cutoffs are fairly relaxed so I don't think we're missing much.
#### Making them more relaxed tended to make the results look weird and not in an interesting way.

## Initial analyses ##
#### MeanVarPlot is depreciated, use `FindVariableGenes`
seurat_object_noBatch <- FindVariableGenes(object = seurat_object_noBatch,
                                     mean.function = ExpMean,
                                     dispersion.function = LogVMR,
                                     x.low.cutoff = 2,
                                     y.cutoff = 0.1)


#### How many genes does this give us?
length(seurat_object_noBatch@var.genes)
#### Less than what was said on the original analysis (4,573)

#### Perform PCA using these genes to see what the populations look like
seurat_object_noBatch <- RunPCA(seurat_object_noBatch, do.print = F)

# Find Outlier cells
PCAPlot(seurat_object_noBatch, 1, 2)
#### There are outliers, which make it difficult to see the strucutre in the main population.
#### These are likely technical artefacts, so we'll remove them.

rotations <- as.data.frame(seurat_object_noBatch@dr$pca@cell.embeddings)
outlier_cells <- rotations[(rotations$PC1<(-20) | rotations$PC2>(15)), ]

tpms_without_outliers <- tpms[,!colnames(tpms) %in% rownames(outlier_cells)]

dim(tpms)
dim(tpms_without_outliers)