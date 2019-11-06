library(qusage)

## Read in the Genesets
All_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Output/Cyto/msigdb.v6.2.entrez.gmt")
KEGG_gmt <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Output/Cyto/c2.cp.kegg.v6.2.entrez.gmt")
GO_terms <- read.gmt("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Output/Cyto/c5.all.v6.2.entrez.gmt")


## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_terms)]
filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_terms) | names(All_gmt) %in% names(KEGG_gmt)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (4326 genesets)
geneset_sizes <- unlist(lapply(filter_gmt, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- filter_gmt[geneset_indices]

CYTOTOX_DATASETS <- filtered_set[grepl("T_CELL_MEDIATED_CYTOTO", names(filtered_set))]
T_CELL_CYTO <- Reduce(intersect, CYTOTOX_DATASETS)

try <- filtered_set[grepl("NATURAL_KILLER_CELL_MEDIATED_IMMUNITY", names(filtered_set))]
NK <- Reduce(intersect, try)

NK


load("./Bulk/Final_edgeR.RData")

library(gplots)
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")[v$targets$group]
data.frame(v$targets$group, col.cell1)
mycol <- colorpanel(1000,"blue","white","red")

# pdf("../Figures/VD2/426_genes_VD2.pdf", width = 10, height = 50)
heatmap.2(v$E[rownames(v$E) %in% T_CELL_CYTO,], scale = "row",
          labRow = ifelse((is.na(v$genes$SYMBOL[T_CELL_CYTO])), v$genes$ENTREZID[T_CELL_CYTO], v$genes$SYMBOL[T_CELL_CYTO]), labCol = NA,
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8, 6), lhei = c(2, 10), 
          hclustfun = hclustAvg, ColSideColors = col.cell1)
dev.off()

v$E[i,]
v$E[T_CELL_CYTO, ]


v$E[rownames(v$E) == "8772"]

v$genes$SYMBOL[T_CELL_CYTO]

v$genes$SYMBOL[rownames(v$E) %in% T_CELL_CYTO]

T_CELL_CYTO <- as.integer(T_CELL_CYTO)
str(i)
str(T_CELL_CYTO)
