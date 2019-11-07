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


## Performing differential gene expression ##
#### Load the .RData object from previous analysis
load("Bulk/Final_edgeR.RData")
## Plotting Most Variable genes ##

#### Get the gene names for the top 100 most variable genes
#### Estimate variance for all genes
#### Remove VD2 from analysis
DF <- as.data.frame(lcpm)
DFa <- DF[, !grepl( "VD2" , names(DF) )]
var_genes <- apply(DFa, 1, var)

#### Remove both TRAV and TRBV genes (expect these to be differentially expressed)
var_genes_df <- as.data.frame(var_genes)
var_genes_df1 <- rownames_to_column(var_genes_df, var = "ENTREZID")
var_genes_df2 <- merge(var_genes_df1, x$genes, by = "ENTREZID")
var_genes_df3 <- var_genes_df2[!grepl("^TRAV", var_genes_df2$SYMBOL),]
var_genes_df4 <- var_genes_df3[!grepl("^TRBV", var_genes_df3$SYMBOL),]
var_genes2 <- extract2(var_genes_df4, 'var_genes') %>% set_names(var_genes_df4$ENTREZID)

#### Select the 150 most variable genes
select_var <- names(sort(var_genes2, decreasing = T))[1:150]

#### Subset logcounts matrix
highly_variable_lcpm <- DFa[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#### Define colour palette
No_VD2 <- droplevels(x$samples$group[x$samples$group != "VD2"])
col.cell1 <- c("#999999","#56B4E9","#E69F00","#009E73")[No_VD2]
mycol <- colorpanel(1000,"blue","white","red")
#### Get the SYMBOLS instead of ENTREZ ID
HV_lpcm <- rownames_to_column(as.data.frame(highly_variable_lcpm), var = "ENTREZID")
HV_lpcm1 <- merge(HV_lpcm, x$genes, by = "ENTREZID")
highly_variable_lcpm_sym <- within(HV_lpcm1, rm(TXCHROM, ENTREZID))

# Plot the heatmap
# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/Heatmap_top150_variable_genes.png",
#     width = 300, height = 350, units = "mm", res = 300)
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
heatmap.2(as.matrix(highly_variable_lcpm_sym[, names(highly_variable_lcpm_sym) != "SYMBOL"]),
          col = mycol,
          trace = "none",
          density.info = "none",
          cex.main = 1.5,
          ColSideColors = col.cell1,
          scale = "none",
          margin = c(10,5), lhei = c(2,10),
          labCol = colnames(highly_variable_lcpm_sym),
          hclustfun = hclustAvg
          #,          labRow = highly_variable_lcpm_sym$SYMBOL
          )
par(xpd = T)
legend(x = 0.87, y = 1.065,
       fill = c("#999999", "#56B4E9",
                "#E69F00", "#009E73",
                "#CC79A7"),
       legend = levels(x$samples$group))
dev.off()


## Differential gene expression analysis ##

#### Create the model and contrast matrix
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

#### Removing heteroscedascity
# biocLite("limma", dependencies = T)
library(limma)
# par(mfrow = c(1,2))
v <- voom(x, design, plot = F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = T)
# plotSA(efit, main = "Final model: Mean−variance trend")
summary(decideTests(efit))

#### Apply a log fold change of 2.5 (More stringently selected DE genes)
tfit <- treat(vfit, lfc = log2(2.5))
dt <- decideTests(efit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,8]!=0)

#### Write the common genes
Common_genes <- tfit$genes$ENTREZ[de.common]

length(Common_genes)
# Entrez <- tfit$genes$ENTREZ[de.common]
# Symbol <- tfit$genes$SYMBOL[de.common]
# Common_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(Common_DE_genes)

de.VD1 <- which(dt[,1]!=0 & dt[,8]==0)
# Entrez <- tfit$genes$ENTREZ[de.VD1]
# Symbol <- tfit$genes$SYMBOL[de.VD1]
# VD1_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(VD1_DE_genes)

de.CD8 <- which(dt[,1]==0 & dt[,8]!=0)
# Entrez <- tfit$genes$ENTREZ[de.CD8]
# Symbol <- tfit$genes$SYMBOL[de.CD8]
# CD8_DE_genes <- cbind(Symbol, Entrez)
# writeCsvO(CD8_DE_genes)

## Up, Down and Both
comm <- length(Common_genes)
VD1.de <- length(de.VD1)
CD8.de <- length(de.CD8)


effector <- which(dt[,5]==1 | dt[,5]==-1)
Entrez <- tfit$genes$ENTREZ[effector]
Symbol <- tfit$genes$SYMBOL[effector]
effector_DE_genes <- cbind(Symbol, Entrez)
writeCsvO(effector_DE_genes)

library(eulerr)
"http://eulerr.co/" #HERE!!

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
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999")[v1$targets$group]
data.frame(v1$targets$group, col.cell1)

pdf("Bulk/Figures/Heatmap_429.pdf", width = 50, height = 10)
heatmap.2(t(v1$E[i,]), scale = "column",
          labCol = ifelse((is.na(v1$genes$SYMBOL[i])), v1$genes$ENTREZID[i], v1$genes$SYMBOL[i]), labRow = NA,
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), hclustfun = hclustAvg, RowSideColors = col.cell1)
dev.off()

# PCA of all shared differentially expressed genes
## Change all 'v' to 'v1' to remove VD2 from PCA
library(ggbiplot)
i <- which(v1$genes$ENTREZID %in% Common_genes)
this <- as.data.frame(v1$E[i,])
that <- rownames_to_column(this, var = "Gene")
head(that)

that1 <- that %>% gather(contains("."), key = "Sample", value = "value")
that2 <- spread(that1, key = Gene, value = value)

# Label groups
subgroup <- ifelse(grepl("VD1.CD27LO", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", that2$Sample), "VD1.CD27LO",
                   ifelse(grepl("CD8.EMRA", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", that2$Sample), "CD8.EMRA",
                          ifelse(grepl("CD8.Naive", that2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", that2$Sample), "CD8.Naive",
                                 ifelse(grepl("VD1.CD27HI", that2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", that2$Sample), "VD1.CD27HI",
                                        ifelse(grepl("VD2", that2$Sample), "VD2", "none")))))
that2$subgroup <- as.factor(subgroup)

# Remove identifiers
that3 <- data.frame(that2[, names(that2) != "Sample"], row.names = that2[, names(that2) == "Sample"])
colnames(that3) <- gsub("^X", "",  colnames(that3))

prin_comp <- prcomp(that3[, names(that3) != "subgroup"], scale. = T)
subgroup <- that3[, "subgroup"]

detach("package:randomForest", unload = T)
dev.off()
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = subgroup, ellipse = F,
              circle = T,
              var.axes = F
) + scale_color_manual(values = cbcols) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
g
g <- g + theme(legend.position = c(0.5, 0.5),
               legend.background = element_rect(size = 0.2,
                                                linetype = "solid", 
                                                colour = "black"),
               legend.text = element_text(size = 6))

g <- g + theme(legend.key.height = unit(.6, "line")) + 
  theme(legend.title = element_blank())

g <- g + guides(colour = guide_legend(override.aes = list(size = 3, linetype = 0, pch = 15)))

ggsave("PCA of shared DE genes (Naive vs Effector).pdf" , plot = g, device = "pdf",
       path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures", scale = .5)


# How many PC's are important?
library(factoextra)
fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue")

# Find the Eigenvalues
eig <- (prin_comp$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

## Store the variances as a dataframe
### Write as a DF
eigenvalues_RNAseq <- data.frame(eig = eig, variance = variance,
                                 cumvariance = cumvar)
# writeCsvO(eigenvalues_RNAseq)

# Store the variances
var <- get_pca_var(prin_comp)

## Find the coordinates of variables
### var$coord[, 1:5]

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation

# ### Find orientation of loadings
loads <- as.data.frame(loadings)
loads[order(loads$PC1, decreasing = T)[1:10],]
loads1 <- tibble:: rownames_to_column(loads, "Parameter")
loads_contrib <- droplevels(subset(loads1, PC1 < -0.26))
loads_contrib$Parameter <- as.factor(loads_contrib$Parameter)

sdev <- prin_comp$sdev

### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
PC1_c <- contrib.var[order(contrib.var$PC1, decreasing = T)[1:50],]

PC1_co <- rownames_to_column(PC1_c, var = "gene")
top100PC1 <- PC1_co[, c("gene")]
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999")[v1$targets$group]
data.frame(v1$targets$group, col.cell1)

pdf("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/PC_Heatmap_top100_contrib_genes_horiz.pdf")
heatmap.2(t(v1$E[top100PC1,]), scale = "column",labRow = FALSE, #v1$genes$SYMBOL[v1$genes$ENTREZID %in% top100PC1],
           labCol = v1$genes$SYMBOL[v1$genes$ENTREZID %in% top100PC1], #colnames(v1), #can also do labCol = groups
          col = mycol, trace = "none", density.info = "none", RowSideColors = col.cell1, key = F,
          margin = c(8,6), lhei = c(2,10))
dev.off()

pdf("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/PC_Heatmap_top100_contrib_genes_vert.pdf")
heatmap.2(v1$E[top100PC1,], scale = "row", labRow = v1$genes$SYMBOL[v1$genes$ENTREZID %in% top100PC1],
          labCol = FALSE, #colnames(v1), #can also do labCol = groups
          col = mycol, trace = "none", density.info = "none", ColSideColors = col.cell1, key = F,
          margin = c(8,6), lhei = c(2,10))
dev.off()

# View comparisons
CD27LO.vs.CD27HI <- topTreat(efit, coef = 1, n = Inf)
# write.csv(CD27LO.vs.CD27HI, file = "./Bulk/Output/CD27LO.vs.CD27HI.csv")

sig_VD1 <- droplevels(subset(CD27LO.vs.CD27HI, adj.P.Val <= 0.05))
sig_VD1$ENTREZID <- as.factor(sig_VD1$ENTREZID)
length(sig_VD1$ENTREZID)

genes_VD1 <- c(levels(sig_VD1$ENTREZID))
# write.table(genes_VD1, file = "Bulk/Output/genes_VD1.tsv", row.names = F)

CD27HI.vs.EMRA <- topTreat(tfit, coef = 2, n = Inf)

CD27HI.vs.Naive <- topTreat(tfit, coef = 3, n = Inf)
CD27HI.vs.VD2 <- topTreat(tfit, coef = 4, n = Inf)

CD27LO.vs.EMRA <- topTreat(tfit, coef = 5, n = Inf)

sig_VD1 <- droplevels(subset(CD27LO.vs.EMRA, adj.P.Val <= 0.05))
head(sig_VD1)

CD27LO.vs.EMRA[grep("TRDV1", CD27LO.vs.EMRA$SYMBOL), ]

head(CD27LO.vs.EMRA)



sig_VD1$ENTREZID <- as.factor(sig_VD1$ENTREZID)



CD27LO.vs.Naive <- topTreat(tfit, coef = 6, n = Inf)
CD27LO.vs.VD2 <- topTreat(tfit, coef = 7, n = Inf)

EMRA.vs.Naive <- topTreat(tfit, coef = 8, n = Inf)
# write.csv(EMRA.vs.Naive, file = "./Bulk/Output/EMRA.vs.Naive.csv")

EMRA.vs.VD2 <- topTreat(tfit, coef = 9, n = Inf)

Naive.vs.VD2 <- topTreat(tfit, coef = 10, n = Inf)


######### AME
sig <- droplevels(subset(CD27LO.vs.CD27HI, adj.P.Val <= 0.01))
sig$Direction <- as.factor(ifelse((sig$logFC < 0), "CD27HI", "CD27LO"))

CD27HI <- droplevels(subset(sig, Direction == "CD27HI"))
genes_27HI <- CD27HI$ENTREZID

CD27LO <- droplevels(subset(sig, Direction == "CD27LO"))
genes_27LO <- CD27LO$ENTREZID

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# CD27HI
CD27HI_genes <- getBM(attributes = c("transcript_start", "transcript_end", "chromosome_name", "ensembl_gene_id"),
                       filters = "entrezgene", values = genes_27HI, mart = mart)


CD27HI_genes$chromosome_name_right <- paste("chr", CD27HI_genes$chromosome_name, sep = "") 
CD27HI_genes$chromosome_name_right <- as.factor(CD27HI_genes$chromosome_name_right)
CD27HI_genes1 <- CD27HI_genes[!grepl("chrCHR", CD27HI_genes$chromosome_name_right),]

CD27HI_ensembl1 <- CD27HI_genes1[!duplicated(CD27HI_genes1$ensembl_gene_id), ]

CD27HI_ensembl1$ensembl_gene_id <- as.factor(CD27HI_ensembl1$ensembl_gene_id)

CD27HI_ensembl2 <- data.frame(ensembl_gene = character(),
                            transcript_start = double(),
                            transcript_end= double(),
                            chromosome_name_right = character(),
                            stringsAsFactors = F)
c <- 1
for(i in levels(CD27HI_ensembl1$ensembl_gene_id)){
  print(i)
  working <- droplevels(subset(CD27HI_ensembl1, ensembl_gene_id == i))
  avg_start <- mean(working$transcript_start)
  avg_end <- mean(working$transcript_end)
  working$chromosome_name_right <- as.factor(working$chromosome_name_right)
  CD27HI_ensembl2[c, "ensembl_gene"] <- i
  CD27HI_ensembl2[c, "transcript_start"] <- avg_start
  CD27HI_ensembl2[c, "transcript_end"] <- avg_end
  CD27HI_ensembl2[c, "chromosome_name_right"] <- as.character(levels(working$chromosome_name_right))
  c <- c + 1
}

CD27HI_ensembl2$transcript_start_ext <- CD27HI_ensembl2$transcript_start - 2000
CD27HI_ensembl2$transcript_end_ext <- CD27HI_ensembl2$transcript_end + 1000
library(GenomicRanges)

gr_27HI <- GRanges(seqnames = Rle(CD27HI_ensembl2$chromosome_name_right),
                   ranges = IRanges(CD27HI_ensembl2$transcript_start_ext,
                                    end = CD27HI_ensembl2$transcript_end_ext),
                   strand = Rle(strand(c(rep("*", length(CD27HI_ensembl2$chromosome_name_right))))),
                   names = CD27HI_ensembl2$ensembl_gene)

df_27HI <- data.frame(seqnames = seqnames(gr_27HI),
                      starts = start(gr_27HI)-1,
                      ends = end(gr_27HI),
                      names = gr_27HI$names
                        )

df_27HI1 <- droplevels(subset(df_27HI, names != "ENSG00000189283"))

subset(df_27HI, starts == 63637707)

write.table(df_27HI1, file = "Bulk/Output/CD27HI_list1.bed", quote = F, sep = "\t", row.names = F, col.names = F)


library(Homo.sapiens)
promoters(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), 2000, -2000)


# CD27LO
CD27LO_genes <- getBM(attributes = c("transcript_start", "transcript_end", "chromosome_name", "ensembl_gene_id"),
                      filters = "entrezgene", values = genes_27LO, mart = mart)


CD27LO_genes$chromosome_name_right <- paste("chr", CD27LO_genes$chromosome_name, sep = "") 
CD27LO_genes$chromosome_name_right <- as.factor(CD27LO_genes$chromosome_name_right)
CD27LO_genes1 <- CD27LO_genes[!grepl("chrCHR", CD27LO_genes$chromosome_name_right),]

CD27LO_ensembl1 <- CD27LO_genes1[!duplicated(CD27LO_genes1$ensembl_gene_id), ]

CD27LO_ensembl1$ensembl_gene_id <- as.factor(CD27LO_ensembl1$ensembl_gene_id)

CD27LO_ensembl2 <- data.frame(ensembl_gene = character(),
                              transcript_start = double(),
                              transcript_end= double(),
                              chromosome_name_right = character(),
                              stringsAsFactors = F)
c <- 1
for(i in levels(CD27LO_ensembl1$ensembl_gene_id)){
  print(i)
  working <- droplevels(subset(CD27LO_ensembl1, ensembl_gene_id == i))
  avg_start <- mean(working$transcript_start)
  avg_end <- mean(working$transcript_end)
  working$chromosome_name_right <- as.factor(working$chromosome_name_right)
  CD27LO_ensembl2[c, "ensembl_gene"] <- i
  CD27LO_ensembl2[c, "transcript_start"] <- avg_start
  CD27LO_ensembl2[c, "transcript_end"] <- avg_end
  CD27LO_ensembl2[c, "chromosome_name_right"] <- as.character(levels(working$chromosome_name_right))
  c <- c + 1
}

CD27LO_ensembl2$transcript_start_ext <- CD27LO_ensembl2$transcript_start - 2000
CD27LO_ensembl2$transcript_end_ext <- CD27LO_ensembl2$transcript_end + 1000
library(GenomicRanges)

gr_27LO <- GRanges(seqnames = Rle(CD27LO_ensembl2$chromosome_name_right),
                   ranges = IRanges(CD27LO_ensembl2$transcript_start_ext,
                                    end = CD27LO_ensembl2$transcript_end_ext),
                   strand = Rle(strand(c(rep("*", length(CD27LO_ensembl2$chromosome_name_right))))),
                   names = CD27LO_ensembl2$ensembl_gene)

df_27LO <- data.frame(seqnames = seqnames(gr_27LO),
                      starts = start(gr_27LO)-1,
                      ends = end(gr_27LO),
                      names = gr_27LO$names
)

write.table(df_27LO, file = "Bulk/Output/CD27LO_list.bed", quote = F, sep = "\t", row.names = F, col.names = F)

### end ####

length(df_27LO$seqnames)

write.table(EMRA.vs.Naive, "Bulk/Output/EMRA_vs_Naive_DGE.txt", sep = "\t", quote = F, row.names = F)

# These are the outputs of EdgeR!!
# head(CD27LO.vs.EMRA)
# head(CD27HI.vs.VD2)
# head(CD27HI.vs.EMRA)
# head(CD27HI.vs.Naive)

## Highlight this for notable genes... (CD28 one example)
pdf("Bulk/Figures/Volcano_Bulk.pdf")
plotMD(efit, column = 1, status = dt[,1], main = colnames(tfit)[1], 
       xlim = c(-2,13), col = c("#009E73", "#999999"))
legend("topright", fill = c("#999999", "black", "#009E73"),
       legend = c("VD1 CD27LO", "No Change", "VD1 CD27HI"))

dev.off()
glMDPlot(tfit, coef = 1, status = dt, main = colnames(efit)[1],
         side.main = "ENTREZID", counts = x$counts, groups = group, launch = T)

# Heatmaps
## Change the top genes, for various heatmaps looking at different sets
library(gplots)

# CD27HI versus CD27LO
CD27LO.vs.CD27HI1 <- CD27LO.vs.CD27HI[!grepl("^TRAV", CD27LO.vs.CD27HI$SYMBOL), ]
CD27LO.vs.CD27HI2 <- CD27LO.vs.CD27HI1[!grepl("^TRBV", CD27LO.vs.CD27HI1$SYMBOL), ]
CD27LO.vs.CD27HI2 <- Take_Sigs(CD27LO.vs.CD27HI2)

CD27LO.vs.CD27HI.topgenes <- CD27LO.vs.CD27HI2$ENTREZID[1:100]
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

# CD27HI versus EMRA
CD27HI.vs.EMRA1 <- CD27HI.vs.EMRA[!grepl("^TRAV", CD27HI.vs.EMRA$SYMBOL),]
CD27HI.vs.EMRA2 <- CD27HI.vs.EMRA1[!grepl("^TRBV", CD27HI.vs.EMRA1$SYMBOL),]
CD27HI.vs.EMRA2 <- Take_Sigs(CD27HI.vs.EMRA2)

CD27HI.vs.EMRA.topgenes <- CD27HI.vs.EMRA2$ENTREZID[1:100]
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

# CD27HI versus Naive
CD27HI.vs.Naive1 <- CD27HI.vs.Naive[!grepl("^TRAV", CD27HI.vs.Naive$SYMBOL),]
CD27HI.vs.Naive2 <- CD27HI.vs.Naive1[!grepl("^TRBV", CD27HI.vs.Naive1$SYMBOL),]
CD27HI.vs.Naive2 <- Take_Sigs(CD27HI.vs.Naive2)

CD27HI.vs.Naive.topgenes <- CD27HI.vs.Naive2$ENTREZID[1:100]
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


# CD27HI versus VD2
CD27HI.vs.VD21 <- CD27HI.vs.VD2[!grepl("^TRAV", CD27HI.vs.VD2$SYMBOL),]
CD27HI.vs.VD22 <- CD27HI.vs.VD21[!grepl("^TRBV", CD27HI.vs.VD21$SYMBOL),]
CD27HI.vs.VD22 <- Take_Sigs(CD27HI.vs.VD22)

CD27HI.vs.VD2.topgenes <- CD27HI.vs.VD22$ENTREZID[1:100]
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


# CD27LO versus VD2
CD27LO.vs.VD21 <- CD27LO.vs.VD2[!grepl("^TRAV", CD27LO.vs.VD2$SYMBOL),]
CD27LO.vs.VD22 <- CD27LO.vs.VD21[!grepl("^TRBV", CD27LO.vs.VD21$SYMBOL),]
CD27LO.vs.VD22 <- Take_Sigs(CD27LO.vs.VD22)

CD27LO.vs.VD2.topgenes <- CD27LO.vs.VD22$ENTREZID[1:100]
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


# CD27LO versus EMRA
CD27LO.vs.EMRA1 <- CD27LO.vs.EMRA[!grepl("^TRAV", CD27LO.vs.EMRA$SYMBOL),]
CD27LO.vs.EMRA2 <- CD27LO.vs.EMRA1[!grepl("^TRBV", CD27LO.vs.EMRA1$SYMBOL),]
CD27LO.vs.EMRA2 <- Take_Sigs(CD27LO.vs.EMRA2)

CD27LO.vs.EMRA.topgenes <- CD27LO.vs.EMRA2$ENTREZID[1:100]
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


# CD27LO versus Naive
CD27LO.vs.Naive1 <- CD27LO.vs.Naive[!grepl("^TRAV", CD27LO.vs.Naive$SYMBOL),]
CD27LO.vs.Naive2 <- CD27LO.vs.Naive1[!grepl("^TRBV", CD27LO.vs.Naive1$SYMBOL),]
CD27LO.vs.Naive2 <- Take_Sigs(CD27LO.vs.Naive2)

CD27LO.vs.Naive.topgenes <- CD27LO.vs.Naive2$ENTREZID[1:100]
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


# EMRA versus Naive
EMRA.vs.Naive1 <- EMRA.vs.Naive[!grepl("^TRAV", EMRA.vs.Naive$SYMBOL),]
EMRA.vs.Naive2 <- EMRA.vs.Naive1[!grepl("^TRBV", EMRA.vs.Naive1$SYMBOL),]
EMRA.vs.Naive2 <- Take_Sigs(EMRA.vs.Naive2)

EMRA.vs.Naive.topgenes <- EMRA.vs.Naive2$ENTREZID[1:100]
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


# EMRA versus VD2
EMRA.vs.VD21 <- EMRA.vs.VD2[!grepl("^TRAV", EMRA.vs.VD2$SYMBOL),]
EMRA.vs.VD22 <- EMRA.vs.VD21[!grepl("^TRBV", EMRA.vs.VD21$SYMBOL),]
EMRA.vs.VD22 <- Take_Sigs(EMRA.vs.VD22)

EMRA.vs.VD2.topgenes <- EMRA.vs.VD22$ENTREZID[1:100]
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


# Naive versus VD2
Naive.vs.VD21 <- Naive.vs.VD2[!grepl("^TRAV", Naive.vs.VD2$SYMBOL),]
Naive.vs.VD22 <- Naive.vs.VD21[!grepl("^TRBV", Naive.vs.VD21$SYMBOL),]
Naive.vs.VD22 <- Take_Sigs(Naive.vs.VD22)

Naive.vs.VD2.topgenes <- Naive.vs.VD22$ENTREZID[1:100]
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


# Shared top 100 DE genes
k <- which(CD27LO.vs.CD27HI.topgenes %in% EMRA.vs.Naive.topgenes)
j <- CD27LO.vs.CD27HI.topgenes[k]
i <- which(v1$genes$ENTREZID %in% j)
mycol <- colorpanel(1000,"blue","white","red")

col.cell1 <- c("#56B4E9","#E69F00","#009E73","#999999")[v1$targets$group]
data.frame(v1$targets$group, col.cell1)

# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/Shared_top100_DEgenes.png",
#     width = 300, height = 300, units = "mm", res = 300)
heatmap.2(v1$E[i,], scale = "row",
          labRow = v1$genes$SYMBOL[i], labCol = group1, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), ColSideColors = col.cell1)
par(xpd = T)
legend(x = 0.87, y = 1.05, 
       fill = c("#999999", "#56B4E9",
                "#E69F00", "#009E73"),
       legend = levels(v1$targets$group))
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
CD27LO.vs.CD27HI.topgenes

# Hallmark_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Hallmark_genesets.rdata")
idx <- ids2indices(Hs.H, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.01))
sig.CD27LO.vs.CD27HI

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
sig.EMRA.vs.Naive




length(Hs.H[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]])


# KEGG
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/kegg_human.rdata")

## Geneset enrichment
idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)
sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.05))
KEGG_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "KEGG_Pathway")
KEGG_GD1 <- KEGG_GD[(KEGG_GD$KEGG_Pathway %in% KEGG_CD8$KEGG_Pathway), ]


head(KEGG_GD)
KEGG_GD$Oneminus <- 1 - KEGG_GD$FDR


p <- ggplot(KEGG_GD, aes(x = KEGG_Pathway, y = Oneminus)) +
  geom_bar(stat = "identity") + coord_flip()
p

kegg_sizes<- as.data.frame(lengths(kegg_human)) %>% rownames_to_column(., var = "kegg_pathway")
colnames(kegg_sizes) <- c("KEGG_Pathway", "Size")
head(kegg_sizes)
GD_K <- merge(KEGG_GD, kegg_sizes, by = "KEGG_Pathway")
GD_K$Num <- GD_K$NGenes/GD_K$Size


Comparison <- "GD"
KEGG_GD2 <- cbind(KEGG_GD, Comparison)
KEGG_GD2

idx <- ids2indices(kegg_human, id = rownames(v))
cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
KEGG_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)


idx <- ids2indices(kegg_human, id = rownames(v))
cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
KEGG_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)


 
rbind(KEGG_GD2, KEGG_CD8a)

KEGG <- rbind(KEGG_CD8a, KEGG_GD2)

KEGG$FDR_Dir <- ifelse((KEGG$Direction == "Down"), KEGG$FDR * -1, KEGG$FDR * 1)

KEGG1 <- KEGG[, c("KEGG_Pathway", "Comparison", "FDR_Dir")]

writeCsvO(KEGG1)

KEGG2 <- spread(KEGG1, key = "KEGG_Pathway", value = "FDR_Dir")

KEGG3 <- KEGG2 %>% remove_rownames %>% column_to_rownames(var = "Comparison")

heatmap.2(t(as.matrix(KEGG3)), col = mycol, trace = "none", density.info = "none", scale = "row"
)
mycol <- colorpanel(1000,"blue","white","red")

mycol



writeCsvO(try)
head(try)

thesePath <- try$KEGG_Pathway








 this <- as.matrix(kegg_human)
idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.EMRA <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsEMRA"])
sig.CD27LO.vs.EMRA <- droplevels(subset(cam.CD27LO.vs.EMRA, FDR <= 0.01))


idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27HI.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27HIvsNaive"])
sig.CD27HI.vs.Naive <- droplevels(subset(cam.CD27HI.vs.Naive, FDR <= 0.01))

head(cam.CD27HI.vs.Naive, 5)

View(contr.matrix)
par(mfrow = c(1,1))
barcodeplot(efit$t[, 1], index = idx$`hsa04650 Natural killer cell mediated cytotoxicity`, main = "Natural killer cell mediated cytotoxicity", labels = c("VD1.CD27HI", "VD1.CD27LO"))




#### Less Interesting Genesets ####
## GO_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/GO_genesets.rdata")

## Remove genesets over or equal to 300 length and below 10 genes
GO_sizes<- as.data.frame(lengths(Hs.c5)) %>% rownames_to_column(., var = "GO_pathway")
colnames(GO_sizes) <- c("GO_Pathway", "Size")
remove_these <- droplevels(subset(GO_sizes, Size >= 300 & Size <= 10))$GO_Pathway
isNameInIndex <- names(Hs.c5) %in% remove_these
GO_terms <- Hs.c5[!isNameInIndex]


idx <- ids2indices(GO_terms, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(camera_results, 5)
sig.CD27LO.vs.CD27HI <- droplevels(subset(camera_results, FDR <= 0.001))
GO_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "GO_Pathway")


GO_sizes <- as.data.frame(lengths(Hs.c5)) %>% rownames_to_column(., var = "GO_pathway")
colnames(GO_sizes) <- c("GO_Pathway", "Size")
head(GO_sizes)
GD_G <- merge(GO_GD, GO_sizes, by = "GO_Pathway")
head(GD_G)

dim(GD_G)

GD_G$Num <- GD_G$NGenes/GD_G$Size
head(GD_G)

GD_G

# ord_GD_GO <- arrange(GD_G, FDR)
# ord_GD_GO$min_log10_FDR <- -log10(ord_GD_GO$FDR)
# 
# head(ord_GD_GO)
# p <- ggplot(ord_GD_GO, aes(x = GO_Pathway, y = min_log10_FDR)) +
#   geom_bar(stat = "identity") 
# 
# p
# ?sec_axis
# Comparison <- "GD"
# KEGG_GD2 <- cbind(KEGG_GD, Comparison)
# KEGG_GD2

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.001))
GO_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
dim(GO_CD8)



Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)





head(GD_GENESETS)

common_elements <- combn(GD_GENESETS, 2, 
                         FUN = function(x) intersect(x[[1]], x[[2]]), simplify = F)

names(common_elements) <- vapply(combn(names(GD_GENESETS), 2, simplify = F), 
                                 paste, collapse = "___", FUN.VALUE = character(1))

head(common_elements)

common_num <- lengths(common_elements)

library(tidyverse)
try <- as.data.frame(common_num)
head(try)
try1 <- rownames_to_column(try, var = "geneset")

this <- try1 %>% separate("geneset", into = c("data1", "data2"), sep = "___")
head(this)

length_of_genesets <- data.frame(data1 = character(),
                                 length1 = double(),
                                 stringsAsFactors = F)
c <- 1
for(i in names(GD_GENESETS)){
  work <- GD_GENESETS[[i]]
  that <- length(work)
  length_of_genesets[c, "data1"] <- i
  length_of_genesets[c, "length1"] <- that
  c <- c + 1
}

head(length_of_genesets)
names(GD_GENESETS) %in% levels(as.factor(this$data1))

this1 <- merge(this, length_of_genesets, by = "data1", all.x = T)
subset(this1, data2 == "GO_RECEPTOR_INHIBITOR_ACTIVITY")


length_of_genesets <- data.frame(data2 = character(),
                                 length2 = double(),
                                 stringsAsFactors = F)
c <- 1
for(i in names(GD_GENESETS)){
  work <- GD_GENESETS[[i]]
  that <- length(work)
  length_of_genesets[c, "data2"] <- i
  length_of_genesets[c, "length2"] <- that
  c <- c + 1
}

this2 <- merge(this1, length_of_genesets, by = "data2", all.x = T)

head(this2)
Geneset_overlap <- this2

Geneset_overlap$specific_1 <- Geneset_overlap$length1 - Geneset_overlap$common_num
Geneset_overlap$specific_2 <- Geneset_overlap$length2 - Geneset_overlap$common_num

Geneset_overlap$neg_agreement <- 2350 - (Geneset_overlap$common_num + Geneset_overlap$specific_1 + Geneset_overlap$specific_2)
head(Geneset_overlap)

Geneset_overlap$kappa <- (Geneset_overlap$common_num + 0)/(Geneset_overlap$common_num + Geneset_overlap$specific_1 +Geneset_overlap$specific_2 + 0)
head(Geneset_overlap)

test <- droplevels(subset(Geneset_overlap, data1 == "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY")) #- NA are where there aren't compared
as.factor(test$data2)

test <- droplevels(subset(Geneset_overlap, data2 == "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY")) #- NA are where there aren't compared
as.factor(test$data1)

writeCsvO(Geneset_overlap)

kappa_scores_GD <- Geneset_overlap[, c("data1", "data2", "kappa")]

library(GGally)
library(network)
mm.net <- network(test[,1:2], directed = T)

ggnet2(mm.net,
       labelon = TRUE,
       size = 2, vjust = -0.6, mode = "kamadakawai", label.size = 3)


test <- droplevels(subset(Geneset_overlap, kappa != 0))

mm.net

install.packages("geomnet")
data(madmen, package = 'geomnet')

head(madmen)

?ggnetworkmap()


KS_GD <- spread(kappa_scores_GD, "data2", value = "kappa")
head(KS_GD)
str(KS_GD)

KS_GD1 <- column_to_rownames(KS_GD, var = "data1")
KS_GD1
is.na(KS_GD1) <- 0

getwd()
GO_GD <- read.csv("Bulk/Output/DEgenes/GO_GD.csv")

head(GO_GD)

GD <- as.character(GO_GD$GO_Geneset)

GD1 <- Geneset_overlap[Geneset_overlap$data1 %in% GD, ]
GD2 <- GD1[GD1$data2 %in% GD, ]

head(GD2)


## Making of a binary gene-term matrix
### Take only the terms that are significant in GD


GD_GENESETS <- Hs.c5[names(Hs.c5) %in% GD]
that <- as.data.frame(unlist(GD_GENESETS)) %>% rownames_to_column(var = "Geneset")
that$Geneset <- gsub("[0-9]{1,4}$", "", that$Geneset)
colnames(that) <- c("Geneset", "Genes")
length(that$Genes)

gene_term_binary <- as.data.frame.matrix(table(that))

head(gene_term_binary)





library(Matrix)
#https://stackoverflow.com/a/51421029/1412059
fun <- function(x) {
  n <- 0.5 + sqrt(0.25 + 2 * length(x)) #calculate dimension
  i <- sequence(seq_len(n - 1)) #row indices
  j <- rep(seq_len(n - 1), seq_len(n - 1)) + 1 # column indices
  sparseMatrix(i = i, j = j, x = x, triangular = TRUE, dims = c(n, n))
}

output <- fun(common_num)
diag(output) <- lengths(GD_GENESETS)
dimnames(output) <- rep(list(names(GD_GENESETS)), 2)

View(output)

head(output)

output1 <- as.data.frame(output)




idx <- ids2indices(Hs.c5, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.01))
GO_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "GO_Geneset")
writeCsvO(GO_GD)

barcodeplot(efit$t[, 1], index = idx$GO_CELL_KILLING, main = "CD27LO vs CD27HI")

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.CD27LO.vs.CD27HI, 5)

sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
GO_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "GO_Geneset")
writeCsvO(GO_CD8)



nrow(GO_GD)
GO_GD$Directed_FDR <- ifelse((GO_GD$Direction == "Up"), GO_GD$FDR * 1, GO_GD$FDR * -1)

head(GO_GD)


ggplot(GO_GD, aes(x = GO_Geneset, y = Directed_FDR)) + geom_point()

# Immunological Signatures
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Immunological_Signatures.rdata")
idx <- ids2indices(Hs.c7, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 10)

barcodeplot(efit$t[, 1], index = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP , 
            index2 = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN, main = "NaiveVsCD8")










install.packages("DOSE")
BiocManager::install("PathwaySplice")
library(DOSE)
library(PathwaySplice)

gene.based.table <- makeGeneTable(Hs.c5)
?makeGeneTable
View(featureBasedData)

res <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),method='Wallenius')

# labeling each node by gene set name
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = FALSE)

# labeling each node by gene set index
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = TRUE)

## Not run: 
# illustrates specification of output file directory
# Enable interactive map and label each node by gene set index
enmap <- enrichmentMap(res,n=10,fixed=FALSE, similarity.threshold=0.3,
                       label.node.by.index = TRUE, output.file.dir=tempdir())

enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = FALSE, output.file.dir=tempdir())
