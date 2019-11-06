# ATACSequencing Data
# Path_to_files <- "/Volumes/2018/noyvertb-cruk-bioinformatics/JackMcMurray/ATACSeq/3_peak_called/"
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

# ARE THESE JUST THE CONSENSUS SEQUENCES?!

library(BiocManager)
# BiocManager::install("Rsubread", dependencies = T)
library(Rsubread)
# source("CreateCounts.R")

## Start of old edgeR pipeline
# use edgeR  START!
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR", dependencies = T)
# library(edgeR)
# files <- list.files(path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/ATAC/Data/", pattern = ".txt$")
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/ATAC/Data/")

# # Replenish "x"
# x <- read.delim(files)
# colnames(x)[1] <- "Merged_Peak_ID"
# head(x)
# colnames(x) <- gsub("_treat_pileup.bdg.bedGraph.avg.over.400.bp", "", colnames(x))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
blklist <- import.bed("Blcklist.bed.gz")

library(ChIPseeker)

peak_files <- list.files(path = "/Volumes/2018/noyvertb-cruk-bioinformatics/JackMcMurray/ATACSeq/3_peak_called", pattern = ".narrowPeak$", full.names = T)

peak_list <- list()
library(DT)
library(ChIPQC)


for(i in 1:length(peak_files)){
  openRegionPeaks <- peak_files[i]
  annot_peaks <- annotatePeak(openRegionPeaks, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
  peak_list[[openRegionPeaks]] <- annot_peaks
  setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/ATAC/Figures/")
  pdf(paste0(openRegionPeaks, ".pdf"))
  plotAnnoPie(annot_peaks)
  dev.off()
  }

for(i in names(peak_list)){
  print(i)
  this_name <- peak_list[[i]]
  pdf(paste0(i, "upset", ".pdf"))
  upsetplot(this_name)
  dev.off()
  }


library(rGREAT)

# seqlevelsStyle(paste("/Volumes/2018/noyvertb-cruk-bioinformatic/JackMcMurray/ATACSeq/3_peak_called/", openRegionPeaks, sep = "")) <- "UCSC"

openRegionPeaks <- "VD1-Eff_peaks.narrowPeak"
great_Job <- submitGreatJob(paste0(Path_to_files, openRegionPeaks), species = "hg19")
availableCategories(great_Job)
great_ResultTable = getEnrichmentTables(great_Job, category = "Pathway Data")
names(great_ResultTable)
great_ResultTable[["PANTHER Pathway"]][1:4, ]





peaks <- dir("/Volumes/2018/noyvertb-cruk-bioinformatics/JackMcMurray/ATACSeq/3_peak_called",
             pattern = "*.narrowPeak$", full.names = T)
peaks_ <- peaks[!grepl("Merged_peaks", peaks)]
myPeaks <- lapply(peaks_, ChIPQC:::GetGRanges, simple = T)

names(myPeaks) <- peaks_
names(myPeaks) <- gsub("/Volumes/2018/noyvertb-cruk-bioinformatics/JackMcMurray/ATACSeq/3_peak_called/", "", names(myPeaks))
names(myPeaks) <- gsub("_peaks.narrowPeak", "", names(myPeaks))

Group <- factor(names(myPeaks))
consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")
library(limma)

setwd("../")
blklist <- import.bed("Data/Blcklist.bed.gz")
View(blklist)

consensusToCount <- consensusToCount[!consensusToCount %over% blklist & !seqnames(consensusToCount) %in% 
                                       "chrM"]

as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(matches("Naive")) %>% 
  vennDiagram(main = "Naive overlapping peaks")
dev.off()
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(matches("Eff|EMRA")) %>% 
  vennDiagram(main = "Effector overlapping peaks")
 dev.off()
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(matches("VD2|CD8.EM")) %>% 
  dplyr::select(-matches("CD8.EMRA")) %>%
  vennDiagram(main = "Naive overlapping peaks")

library(tidyr)

data <- as.data.frame(elementMetadata(consensusToCount)) %>% column_to_rownames(var = "consensusIDs") %>% dplyr::select(-CD8.CM, -CD8.EM, -VD2) 
data1 <- as.data.frame(elementMetadata(consensusToCount)) %>% column_to_rownames(var = "consensusIDs") %>% dplyr::select(-CD8.CM, -CD8.EM, -VD2) %>% 
  as.matrix %>% t %>% prcomp

## ROWNAMES ARE NOT THE PEAKS?! CONSENSUS ID = PEAKS. NEED TO MAINTAIN THE CONSENSUS IDS - put them as rownames? CONSENSUS still loses PRDM1?!

data1a <- data1$x %>% data.frame %>% mutate(Group = rownames(.))

data1a$Group <- gsub("VD1.Eff", "VD1.CD27LO", data1a$Group)
data1a$Group <- gsub("VD1.Naive", "VD1.CD27HI", data1a$Group)
data1a$Group <- gsub("CD8.NAIVE", "CD8.Naive", data1a$Group)

library(ggbiplot)
ggbiplot(data1, obs.scale = 1, var.scale = 1, ellipse = T,
              circle = F,
              var.axes = F,
              choices = c(1,2)) + 
  scale_color_manual(values = cbcols) + theme_bw() + theme(legend.direction = "horizontal", 
               legend.position = "top")

ggplot(data1a, aes(x = PC1, y = PC2, colour = Group)) + 
  geom_point(size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + scale_color_manual(values = cbcols) +
  xlab("PC1 (46.3% explained var.)") + ylab("PC2 (33.1% explained var.")

library(factoextra)
fviz_screeplot(data1, ncp = 10, choice = "variance")
eig <- (data1$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

Eigenvalues.All <- data.frame(eig = eig, variance = variance,
                              cumvariance = cumvar)

# Store the variances
var <- get_pca_var(data1)

## Find the correlation between variables and principal components
loadings <- data1$rotation
loads <- as.data.frame(loadings)
load1 <- tibble:: rownames_to_column(loads, "Parameter")
load1$Parameter <- as.factor(load1$Parameter)

# Biggest contributions to a given PC
# load1[order(load1$PC1, decreasing = T)[1:4],]
# droplevels(subset(load1, Parameter == "CD8.10.Class.2.Epithelium"))
# levels(load1$Parameter)
sdev <- data1$sdev

### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
Contrib <- contrib.var %>% rownames_to_column(var = "Peak.ID")

library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
  rowSums

table(occurrences) %>% rev %>% cumsum

### Can't remove those from artefact or those rightly being expressed in one sample due to lack of bio repeat

#

bamsToCount <- dir("/Volumes/2018/noyvertb-cruk-bioinformatics/JackMcMurray/ATACSeq/2_bam/picard_processed", full.names = T, pattern = "*.\\.bam$")

# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = T, 
                           countMultiMappingReads = F, maxFragLength = 400)
myCounts <- fcResults$counts

colnames(myCounts) <- gsub("X.Volumes.2018.noyvertb.cruk.bioinformatics.JackMcMurray.ATACSeq.2_bam.picard_processed.", "", colnames(myCounts))
colnames(myCounts) <- gsub(".bam", "", colnames(myCounts))
colnames(myCounts) <- gsub("VD1.Eff", "VD1.CD27LO", colnames(myCounts))
colnames(myCounts) <- gsub("VD1.Naive", "VD1.CD27HI", colnames(myCounts))
colnames(myCounts) <- gsub("CD8.NAIVE", "CD8.Naive", colnames(myCounts))
library(BSgenome.Hsapiens.UCSC.hg19)

head(myCounts)

str(myCounts)

results(myCounts, c("Group", "VD1.CD27LO", "VD1.CD27HI"), format = "GRanges")

GRa

###
x <- x[, c(1, 20:26)]
x <- column_to_rownames(x, var = "Merged_Peak_ID")
x <- DGEList(counts = x)

samplenames <- substring(colnames(x), nchar(colnames(x)))
samplenames <- gsub("_treat_pileup.bdg.bedGraph.avg.over.400.bp", "", colnames(x))
colnames(x) <- samplenames

# Need to label based on cell type and lane
## Cell types
group <- colnames(x)
x$samples$group <- as.factor(group)


library(rtracklayer)

head(x$samples)

blklist <- import.bed("Blcklist.bed.gz")
openRegionPeaks <- "Merged_peaks.narrowPeak"


MacsCalls_chr20_filteredAnno <- annotatePeak(opd, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)






