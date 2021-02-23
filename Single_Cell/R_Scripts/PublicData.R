# A script to investigate the single-cell RNA sequencing data from Gutierrez study
# Author: Jack McMurray
set.seed(2021) # Year of analysis

# Load up packages and make common variables
cbcols <- c("VD1.CD27LO" = "#999999",
            "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00",
            "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")



SC <- read.delim("./Single_Cell/Data/Gutierrez/GSE124731_single_cell_rnaseq_gene_counts.txt")
SC_meta <- read.delim("./Single_Cell/Data/Gutierrez/GSE124731_single_cell_rnaseq_meta_data.txt")



head(SC)[, 1:10]
head(SC_meta)

vd1 <- droplevels(subset(SC_meta, cell.type == "Vd1"))

nrow(vd1)
