# A script to create counts from GrCh38 aligned bulk RNASequencing data from GammaDelta and CD8 populations
## Author: Jack McMurray
set.seed(2021) # Year of reanalysis

## Packages
library(tidyverse)

# Comparisons and Colours
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"), c("VD1.CD27HI", "CD8.EMRA"), c("VD1.CD27HI", "CD8.Naive"), c("VD1.CD27HI", "VD2"),
                       c("VD1.CD27LO", "CD8.EMRA"), c("VD1.CD27LO", "CD8.Naive"), c("VD1.CD27LO", "VD2"),
                       c("CD8.EMRA", "CD8.Naive"), c("CD8.EMRA", "VD2"), c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999", "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00", "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")

# BiocManager::install("Rsubread", dependencies = T)
library(Rsubread)

# Run after having run the terminal protocol
# bam.files <- get_sorbam("/Volumes/ResearchData/Willcox Group/Jack/GD_RNA_Comb/BAM/")
bam.files <- list.files("/Volumes/noyvertb-cruk-bioinformatics/JackMcMurray/GD_RNASeq/bam_38", pattern = ".bam$", full.names = T)

## Summary of the proportion of reads that are mapped
props <- propmapped(files = bam.files)
# props

## Counting
### Contains inbuilt annotation for hg19 genome assembly
fc <- featureCounts(bam.files, annot.inbuilt = "hg38", isPairedEnd = T)
Counts <- as.data.frame(fc$counts)

# Gain genelengths
## Create a dataframe to allow merging
Counts1 <- tibble:: rownames_to_column(Counts, var = "GeneID")
Counts1$GeneID <- as.factor(Counts1$GeneID)
annotation <- as.data.frame(fc$annotation)
ann <- annotation[ , which(names(annotation) %in% c("GeneID", "Length"))]
ann$GeneID <- as.factor(ann$GeneID)

Counts2 <- droplevels(merge(Counts1, ann, by = "GeneID"))

Counts2$Length <- NULL

write.table("./Bulk/Counts/Counts_GD.txt", x = Counts2, row.names = F, quote = F, sep = "\t")

