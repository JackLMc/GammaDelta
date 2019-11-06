# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsubread", dependencies = T)
library(Rsubread)
library(edgeR)
library(Homo.sapiens)
library(Glimma)
source("Functions.R")
# source("CreateCounts.R")

# use edgeR  START!
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR", dependencies = T)
library(edgeR)
files <- list.files(path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Counts", pattern = ".txt$")
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Counts")

# Replenish "x"
x <- readDGE(files, columns = c(1,3))

samplenames <- substring(colnames(x), length(files), nchar(colnames(x)))
colnames(x) <- samplenames

# Need to label based on cell type and lane
## Cell types
group <- ifelse(grepl("VD1.CD27LO", colnames(x)) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", colnames(x)), "VD1.CD27LO",
                ifelse(grepl("CD8.EMRA", colnames(x)) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", colnames(x)), "CD8.EMRA",
                       ifelse(grepl("CD8.Naive", colnames(x)) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", colnames(x)), "CD8.Naive",
                              ifelse(grepl("VD1.CD27HI", colnames(x)) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", colnames(x)), "VD1.CD27HI",
                                     ifelse(grepl("VD2", colnames(x)), "VD2", "none")))))
x$samples$group <- as.factor(group)

## Lane
lane <- ifelse(grepl("L001", colnames(x)), "L001",
               ifelse(grepl("L002", colnames(x)), "L002",
                      ifelse(grepl("L003", colnames(x)), "L003",
                             ifelse(grepl("L004", colnames(x)), "L004", "none"))))
x$samples$lane <- as.factor(lane)

# Rename columns - annoying
colnames(x) <- c("28MD.VD1.CD27LO.L001", "28MD.VD1.CD27LO.L002", "28MD.VD1.CD27LO.L003", "28MD.VD1.CD27LO.L004",
                 "28MD.CD8.EMRA.L001", "28MD.CD8.EMRA.L002","28MD.CD8.EMRA.L003", "28MD.CD8.EMRA.L004",
                 "28MD.CD8.Naive.L001", "28MD.CD8.Naive.L002", "28MD.CD8.Naive.L003", "28MD.CD8.Naive.L004",
                 "28MD.VD1.CD27HI.L001", "28MD.VD1.CD27HI.L002", "28MD.VD1.CD27HI.L003", "28MD.VD1.CD27HI.L004",
                 "28MD.VD2.L001", "28MD.VD2.L002", "28MD.VD2.L003", "28MD.VD2.L004",
                 "31EN.CD8.Naive.L001", "31EN.CD8.Naive.L002", "31EN.CD8.Naive.L003", "31EN.CD8.Naive.L004",
                 "31EN.CD8.EMRA.L001", "31EN.CD8.EMRA.L002", "31EN.CD8.EMRA.L003", "31EN.CD8.EMRA.L004",
                 "31EN.VD1.CD27LO.L001", "31EN.VD1.CD27LO.L002", "31EN.VD1.CD27LO.L003", "31EN.VD1.CD27LO.L004",
                 "31EN.VD1.CD27HI.L001", "31EN.VD1.CD27HI.L002", "31EN.VD1.CD27HI.L003", "31EN.VD1.CD27HI.L004",
                 "31EN.VD2.L001", "31EN.VD2.L002", "31EN.VD2.L003", "31EN.VD2.L004")


# Gaining second dataframe (Symbols)
# biocLite("Homo.sapiens", dependencies = T)
library(Homo.sapiens)
geneid <- rownames(x)
genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM"), 
                keytype = "ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes



# Calculate counts per million, log counts per million (can also do RPKM [function rpkm()])
cpm <- cpm(x)
lcpm <- cpm(x, log = T)

# Remove lowly expressed transcripts
table(rowSums(x$counts==0)==40) # Show me the amount of transcripts that are zero for all 40 samples
keep.exprs <- rowSums(cpm>1)>=4 # Genes must have a cpm above 1 and be expressed in at least 4 runs to be kept
x <- x[keep.exprs,, keep.lib.sizes = F]
dim(x)

## Showing the removal of the lowly expressed transcripts
# library(RColorBrewer)
# nsamples <- ncol(x)
# col <- brewer.pal(nsamples, "Paired")
# par(mfrow = c(1,2))
# plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2, 
#      main = "", xlab = "")
# title(main = "A. Raw data", xlab = "Log-cpm")
# abline(v = 0, lty = 3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col = col[i], lwd = 2)
# }
# legend("topright", samplenames, text.col = col, bty = "n")
# lcpm <- cpm(x, log = T)
# plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2, 
#      main = "", xlab = "")
# title(main = "B. Filtered data", xlab = "Log-cpm")
# abline(v = 0, lty = 3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col = col[i], lwd = 2)
# }
# legend("topright", samplenames, text.col = col, bty = "n")


# Normalising gene expression
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# save(x, file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/x.RData")

## Showing theoretical effect of normalisation
# x2 <- x
# x2$samples$norm.factors <- 1
# x2$counts[,1] <- ceiling(x2$counts[,1]*0.05) # reduce first sample to 5%
# x2$counts[,2] <- x2$counts[,2]*5 # Inflate second sample by x5

### Graph
# par(mfrow = c(1,2))
# lcpm <- cpm(x2, log = T)
# boxplot(lcpm, las = 2, col = col, main = "")
# title(main = "A. Example: Unnormalised data", ylab = "Log-cpm")
# x2 <- calcNormFactors(x2)  
# x2$samples$norm.factors
# lcpm <- cpm(x2, log = T)
# boxplot(lcpm, las = 2, col = col, main = "")
# title(main = "B. Example: Normalised data", ylab = "Log-cpm")

# Unsupervised clustering of samples & runs
# lcpm <- cpm(x, log = T)
# par(mfrow = c(1,2))
# col.group <- as.factor(group)
# levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") 
# col.group <- as.character(col.group)
# col.lane <- as.factor(lane)
# levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
# col.lane <- as.character(col.lane)
# plotMDS(lcpm, labels = group, col = col.group)
# title(main = "A. Sample groups")
# plotMDS(lcpm, labels = lane, col = col.lane, dim = c(1,2))
# title(main = "B. Sequencing lanes")

## Online launch of this.
# biocLite("Glimma", dependencies = T)
library(Glimma)
# glMDSPlot(lcpm, labels = paste(group, lane, sep = "_"), 
#           groups = x$samples[,c(2,5)], launch = F)

# Differential gene expression (START FROM HERE)
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/x.RData")
design <- model.matrix(~0 + group + lane)
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

## Removing heteroscedascity
par(mfrow = c(1,2))
v <- voom(x, design, plot = F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

# Extra Filter
## VD1
x1 <- x
x1$samples <- droplevels(subset(x1$samples, group == "VD1.CD27LO" | group == "VD1.CD27HI"))
x1$counts <- as.data.frame(x1$counts)
x1$counts <- x1$count %>% dplyr:: select(contains("VD1"))
x1$counts <- as.matrix(x1$counts)
group1 <- x1$samples$group
lane1 <- x1$samples$lane
design1 <- model.matrix(~0 + group1 + lane1)
colnames(design1) <- gsub("group1", "", colnames(design1))

cpm1 <- cpm(x1)
lcpm1 <- cpm(x1, log = T)
table(rowSums(x1$counts==0)==40) # Show me the amount of transcripts that are zero for all 40 samples
keep.exprs <- rowSums(cpm1>1)>=4 # Genes must have a cpm above 1 and be expressed in at least 4 runs to be kept
x1 <- x1[keep.exprs,, keep.lib.sizes = F]
v1 <- voom(x1, design1, plot = F)
contr.matrix_VD1 <- makeContrasts(CD27LOvsCD27HI = VD1.CD27LO - VD1.CD27HI,
                                  levels = colnames(design1))

vfit1 <- lmFit(v1, design1)
vfit1 <- contrasts.fit(vfit1, contrasts = contr.matrix_VD1)
efit1 <- eBayes(vfit1)
summary(decideTests(efit1))


tfit1 <- treat(vfit1, lfc = 2)
dt1 <- decideTests(tfit1)
summary(dt1)

## CD8
x2 <- x
View(x2$samples$group)
x2$samples <- droplevels(subset(x2$samples, group == "CD8.EMRA" | group == "CD8.Naive"))
x2$counts <- as.data.frame(x2$counts)
x2$counts <- x2$count %>% dplyr:: select(contains("CD8"))
x2$counts <- as.matrix(x2$counts)
group2 <- x2$samples$group
lane2 <- x2$samples$lane
design2 <- model.matrix(~0 + group2 + lane2)
colnames(design2) <- gsub("group2", "", colnames(design2))

cpm2 <- cpm(x2)
lcpm2 <- cpm(x2, log = T)
table(rowSums(x2$counts==0)==40) # Show me the amount of transcripts that are zero for all 40 samples
keep.exprs <- rowSums(cpm2>1)>=4 # Genes must have a cpm above 1 and be expressed in at least 4 runs to be kept
x2 <- x2[keep.exprs,, keep.lib.sizes = F]
v2 <- voom(x2, design2, plot = F)
contr.matrix_CD8 <- makeContrasts(EMRAvsNaive = CD8.EMRA - CD8.Naive,
                                  levels = colnames(design2))

vfit2 <- lmFit(v2, design2)
vfit2 <- contrasts.fit(vfit2, contrasts = contr.matrix_CD8)
efit2 <- eBayes(vfit2)
summary(decideTests(efit2))


tfit2 <- treat(vfit2, lfc = 2)
dt2 <- decideTests(tfit2)
summary(dt2)

CD8_comp <- as.data.frame(dt2)
VD1_comp <- as.data.frame(dt1)

head(CD8_comp)

CD8a <- rownames_to_column(CD8_comp, var = "Entrez")
VD1a <- rownames_to_column(VD1_comp, var = "Entrez")

head(CD8a)

CD8_filter <- droplevels(subset(CD8a, EMRAvsNaive != 0 ))
VD1_filter <- droplevels(subset(VD1a, CD27LOvsCD27HI != 0 ))


VD1_up <- droplevels(subset(VD1_filter, CD27LOvsCD27HI == 1))
CD8_up <- droplevels(subset(CD8_filter, EMRAvsNaive == 1))
length(CD8_up$Entrez)
merge_up <- merge(VD1_up, CD8_up, by = "Entrez")
length(merge_up$Entrez)
 
shared <- VD1_up$Entrez[VD1_up$Entrez %in% merge_up$Entrez]

VD1_specific_up <- VD1_up$Entrez[!VD1_up$Entrez %in% shared]
CD8_specific_up <- CD8_up$Entrez[!CD8_up$Entrez %in% shared]


Specific_VD1_up_genes  <- listofgenes[listofgenes$ENTREZID %in% VD1_specific_up, ]
Shared_up_genes <- listofgenes[listofgenes$ENTREZID %in% shared, ]

Specific_CD8_up_genes  <- listofgenes[listofgenes$ENTREZID %in% CD8_specific_up, ]
length(Specific_CD8_up_genes$ENTREZID)
writeCsvO(Specific_CD8_up_genes)
writeCsvO(Specific_VD1_up_genes)
writeCsvO(Shared_up_genes)
head(these)
listofgenes <- x$genes



length(CD8_up$Entrez[CD8_up$Entrez %in% merge_up$Entrez])

# Down
VD1_down <- droplevels(subset(VD1_filter, CD27LOvsCD27HI == -1))
CD8_down <- droplevels(subset(CD8_filter, EMRAvsNaive == -1))
length(VD1_down$Entrez)
merge_down <- merge(VD1_down, CD8_down, by = "Entrez")
length(merge_down$Entrez)

length(VD1_up$Entrez[VD1_up$Entrez %in% merge_up$Entrez])
length(CD8_up$Entrez[CD8_up$Entrez %in% merge_up$Entrez])


de.common_filter <- which(dt1[,1]!=0 & dt2[,1]!=0)


# plotSA(efit, main = "Final model: Mean−variance trend")

## Number of differentially expressed genes
summary(decideTests(efit))

### More stringently selected DE genes
tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,8]!=0)

#### Write the common genes
Common_genes <- tfit$genes$ENTREZ[de.common]
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


par(mfrow = c(1,1))
vennDiagram(dt[,c(1,8)], circle.col=c("turquoise", "salmon"))

# Heatmap of shared differentially expressed genes 
## Remove VD2
v1 <- v
check <- as.data.frame(v$E)
colnames(check)
v1$targets <- droplevels(subset(v1$targets, group != "VD2"))
v1$E <- as.data.frame(v1$E)
v1$E <- v1$E %>% dplyr::select(-contains("VD2"))
colnames(v1$E)
v1$E <- as.matrix(v1$E)
j <- !grepl("VD2", group)
group1 <- group[j]
i <- which(v1$genes$ENTREZID %in% Common_genes)


mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v1$E[i,], scale = "row",
          labRow = v1$genes$SYMBOL[i], labCol = colnames(v1), #can also do labCol = groups
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
          main = "Shared DE genes (Naive v Effectors)")
dev.off()

# PCA of all shared differentially expressed genes
## Change all 'v' to 'v1' to remove VD2 from PCA
library(ggbiplot)
i <- which(v$genes$ENTREZID %in% Common_genes)
this <- as.data.frame(v$E[i,])
that <- rownames_to_column(this, var = "Gene")
that1 <- that %>% gather(contains("L00"), key = "Sample", value = "value")
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


g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = subgroup, ellipse = T,
              circle = T,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

g <- g + ggtitle("PCA of shared DE genes (Naive vs Effector)")

ggsave("PCA of shared DE genes (Naive vs Effector).png" , plot = g, device = "png",
       path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)

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
contrib.var[order(contrib.var$PC1, decreasing = T)[1:10],]

# Use the analysis to compare the values for the groups across PC
# Store the contribution to each PC
PC_TCGA <- as.data.frame(prin_comp$x)
PC_TCGA1 <- tibble:: rownames_to_column(PC_TCGA, "Sample")
PC_TCGA1$Sample <- as.factor(PC_TCGA1$Sample)

## Gather the PCs 
PC_TCGA2 <- PC_TCGA1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
PC_TCGA2$Component <- as.factor(PC_TCGA2$Component)

# Plot
## Test for normality
subgroup <- ifelse(grepl("VD1.CD27LO", PC_TCGA2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", PC_TCGA2$Sample), "VD1.CD27LO",
                   ifelse(grepl("CD8.EMRA", PC_TCGA2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", PC_TCGA2$Sample), "CD8.EMRA",
                          ifelse(grepl("CD8.Naive", PC_TCGA2$Sample) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", PC_TCGA2$Sample), "CD8.Naive",
                                 ifelse(grepl("VD1.CD27HI", PC_TCGA2$Sample) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", PC_TCGA2$Sample), "VD1.CD27HI",
                                        ifelse(grepl("VD2", PC_TCGA2$Sample), "VD2", "none")))))
PC_TCGA2$subgroup <- as.factor(subgroup)
Normality <- PC_TCGA2
Normality$uniq <- as.factor(paste(Normality$subgroup, Normality$Component, sep = ","))

## Shapiro test
normal_list <- list()
c <- 1
for(i in levels(Normality$uniq)){
  name <- basename(i)
  cat('Processing', i, '\n')
  working <- droplevels(subset(Normality, uniq == i))
  stat <- shapiro.test(working[,"ComponentScore"])
  normal <- as.data.frame(stat$p.value)
  normal_list[[i]] <- cbind(i, normal)
  c <- c + 1
}
z <- do.call(rbind, normal_list)
rownames(z) <- c()
head(z)

## QQ plots
for(i in levels(Normality$uniq)){
  name <- basename(i)
  cat('Processing', i, '\n')
  working <- droplevels(subset(Normality, uniq == i))
  temp_plot <- ggplot(working) +
    stat_qq(aes(sample = ComponentScore))
  filen <- paste0(i, ".png")
  ggsave(filen, plot = temp_plot, device = "tiff",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Figures/PCA/Normality_PC",
         height=5, width=5, units='in', dpi=600)
}

# Compare across components 
for(i in levels(PC_TCGA2$Component)){
  name <- basename(i)
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(PC_TCGA2, Component == i))
  Chosen$Rank <- rank(Chosen$ComponentScore)
  temp_plot <- ggComponent(Chosen)
  filen <- paste0(i, ".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Figures/PCA/Compare_PC",
         height=6, width=7, units='in', dpi=600)
}

# Combine stats into a list
stat_list <- list()
c <- 1
for(i in levels(PC_TCGA2$Component)) {
  name <- basename(i)
  cat('Processing', i, '\n')
  workingon <- droplevels(subset(PC_TCGA2, Component == i))
  workingon$Rank <- rank(workingon$ComponentScore)
  # assign your ggplot call to the i'th position in the list
  x <- compare_means(Rank ~ subgroup, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}

# Bind and remove row names
z <- do.call(rbind, stat_list)
head(z)
rownames(z) <- c()
RNAseq_PCAStats <- z

# Write out the statistics
# writeCsvO(RNAseq_PCAStats)

# write.fit(tfit, dt, file="results.txt")


# View comparisons
CD27LO.vs.CD27HI <- topTreat(tfit, coef = 1, n = Inf)
CD27HI.vs.EMRA <- topTreat(tfit, coef = 2, n = Inf)
CD27HI.vs.Naive <- topTreat(tfit, coef = 3, n = Inf)
CD27HI.vs.VD2 <- topTreat(tfit, coef = 4, n = Inf)

CD27LO.vs.EMRA <- topTreat(tfit, coef = 5, n = Inf)
CD27LO.vs.Naive <- topTreat(tfit, coef = 6, n = Inf)
CD27LO.vs.VD2 <- topTreat(tfit, coef = 7, n = Inf)

EMRA.vs.Naive <- topTreat(tfit, coef = 8, n = Inf)
EMRA.vs.VD2 <- topTreat(tfit, coef = 9, n = Inf)

Naive.vs.VD2 <- topTreat(tfit, coef = 10, n = Inf)

# These are the outputs of EdgeR!!
# head(CD27LO.vs.CD27HI)
# head(CD27HI.vs.VD2)
# head(CD27HI.vs.EMRA)
# head(CD27HI.vs.Naive)

plotMD(tfit, column = 1, status = dt[,1], main = colnames(tfit)[1], 
       xlim = c(-8,13))

glMDPlot(tfit, coef = 1, status = dt, main = colnames(tfit)[1],
         side.main = "ENTREZID", counts = x$counts, groups = group, launch = F)


# Heatmaps
## Change the top genes, for various heatmaps looking at different sets
library(gplots)

# CD27HI versus CD27LO
CD27LO.vs.CD27HI1 <- CD27LO.vs.CD27HI[!grepl("^TRAV", CD27LO.vs.CD27HI$SYMBOL),]
CD27LO.vs.CD27HI2 <- CD27LO.vs.CD27HI1[!grepl("^TRBV", CD27LO.vs.CD27HI1$SYMBOL),]

CD27LO.vs.CD27HI.topgenes <- CD27LO.vs.CD27HI2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27LO.vs.CD27HI.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
          main = "CD27HI vs CD27LO")

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
CD27HI.vs.EMRA.topgenes <- CD27HI.vs.EMRA2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27HI.vs.EMRA.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27HI vs EMRA")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27HI")
# 
# CD27HI_EMRA_DEgenes <- this
# writeCsvO(CD27HI_EMRA_DEgenes)

# CD27HI versus Naive
CD27HI.vs.Naive1 <- CD27HI.vs.Naive[!grepl("^TRAV", CD27HI.vs.Naive$SYMBOL),]
CD27HI.vs.Naive2 <- CD27HI.vs.Naive1[!grepl("^TRBV", CD27HI.vs.Naive1$SYMBOL),]
CD27HI.vs.Naive.topgenes <- CD27HI.vs.Naive2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27HI.vs.Naive.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27HI vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27HI")
# head(this)
# CD27HI_Naive_DEgenes <- this
# writeCsvO(CD27HI_Naive_DEgenes)


# CD27HI versus VD2
CD27HI.vs.VD21 <- CD27HI.vs.VD2[!grepl("^TRAV", CD27HI.vs.VD2$SYMBOL),]
CD27HI.vs.VD22 <- CD27HI.vs.VD21[!grepl("^TRBV", CD27HI.vs.VD21$SYMBOL),]
CD27HI.vs.VD2.topgenes <- CD27HI.vs.VD22$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27HI.vs.VD2.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27HI vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27HI")
# head(this)
# CD27HI_VD2_DEgenes <- this
# writeCsvO(CD27HI_VD2_DEgenes)


# CD27LO versus VD2
CD27LO.vs.VD21 <- CD27LO.vs.VD2[!grepl("^TRAV", CD27LO.vs.VD2$SYMBOL),]
CD27LO.vs.VD22 <- CD27LO.vs.VD21[!grepl("^TRBV", CD27LO.vs.VD21$SYMBOL),]
CD27LO.vs.VD2.topgenes <- CD27LO.vs.VD22$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27LO.vs.VD2.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27LO vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27LO")
# head(this)
# CD27LO_VD2_DEgenes <- this
# writeCsvO(CD27LO_VD2_DEgenes)


# CD27LO versus EMRA
CD27LO.vs.EMRA1 <- CD27LO.vs.EMRA[!grepl("^TRAV", CD27LO.vs.EMRA$SYMBOL),]
CD27LO.vs.EMRA2 <- CD27LO.vs.EMRA1[!grepl("^TRBV", CD27LO.vs.EMRA1$SYMBOL),]
CD27LO.vs.EMRA.topgenes <- CD27LO.vs.EMRA2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27LO.vs.EMRA.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27LO vs EMRA")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27LO")
# head(this)
# CD27LO_EMRA_DEgenes <- this
# writeCsvO(CD27LO_EMRA_DEgenes)


# CD27LO versus Naive
CD27LO.vs.Naive1 <- CD27LO.vs.Naive[!grepl("^TRAV", CD27LO.vs.Naive$SYMBOL),]
CD27LO.vs.Naive2 <- CD27LO.vs.Naive1[!grepl("^TRBV", CD27LO.vs.Naive1$SYMBOL),]
CD27LO.vs.Naive.topgenes <- CD27LO.vs.Naive2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% CD27LO.vs.Naive.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "CD27LO vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27LO")
# head(this)
# CD27LO_Naive_DEgenes <- this
# writeCsvO(CD27LO_Naive_DEgenes)


# EMRA versus Naive
EMRA.vs.Naive1 <- EMRA.vs.Naive[!grepl("^TRAV", EMRA.vs.Naive$SYMBOL),]
EMRA.vs.Naive2 <- EMRA.vs.Naive1[!grepl("^TRBV", EMRA.vs.Naive1$SYMBOL),]
EMRA.vs.Naive.topgenes <- EMRA.vs.Naive2$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% EMRA.vs.Naive.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "EMRA vs Naive")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "EMRA")
# head(this)
# EMRA_Naive_DEgenes <- this
# writeCsvO(EMRA_Naive_DEgenes)


# EMRA versus VD2
EMRA.vs.VD21 <- EMRA.vs.VD2[!grepl("^TRAV", EMRA.vs.VD2$SYMBOL),]
EMRA.vs.VD22 <- EMRA.vs.VD21[!grepl("^TRBV", EMRA.vs.VD21$SYMBOL),]
EMRA.vs.VD2.topgenes <- EMRA.vs.VD22$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% EMRA.vs.VD2.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "EMRA vs VD2")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "EMRA")
# tail(this)
# EMRA_VD2_DEgenes <- this
# writeCsvO(EMRA_VD2_DEgenes)


# Naive versus VD2
Naive.vs.VD21 <- Naive.vs.VD2[!grepl("^TRAV", Naive.vs.VD2$SYMBOL),]
Naive.vs.VD22 <- Naive.vs.VD21[!grepl("^TRBV", Naive.vs.VD21$SYMBOL),]
Naive.vs.VD2.topgenes <- Naive.vs.VD22$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% Naive.vs.VD2.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), dendrogram = "column", main = "Naive vs VD2")

## gaining a dataframe with the differentially expressed
# this <- Naive.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "Naive")
# head(this)
# Naive_VD2_DEgenes <- this
# writeCsvO(Naive_VD2_DEgenes)


# Shared top 100 DE genes
k <- which(CD27LO.vs.CD27HI.topgenes %in% EMRA.vs.Naive.topgenes)
j <- CD27LO.vs.CD27HI.topgenes[k]
i <- which(v1$genes$ENTREZID %in% j)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v1$E[i,], scale = "row",
          labRow = v1$genes$SYMBOL[i], labCol = group1, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
          main = "Shared top 100 DE genes (Naive v Effectors)")

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








# Geneset Enrichment Analysis
## GO_genesets
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/RData_Objects/GO_genesets.rdata")
idx <- ids2indices(Hs.c5, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 5)
barcodeplot(efit$t[, 1], index = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP , 
            index2 = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN, main = "NaiveVsCD8")

# Immunological Signatures
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/RData_Objects/Immunological_Signatures.rdata")
idx <- ids2indices(Hs.c7, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 10)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, PValue < 0.05))


barcodeplot(efit$t[, 1], index = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP , 
            index2 = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN, main = "NaiveVsCD8")

# Hallmark_genesets
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/RData_Objects/Hallmark_genesets.rdata")
idx <- ids2indices(Hs.H, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, PValue < 0.05))
