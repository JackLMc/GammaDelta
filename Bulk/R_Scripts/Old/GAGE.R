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




# Network analysis
## source("https://bioconductor.org/biocLite.R")
## biocLite("gage", dependencies = T)
## biocLite("pathview", dependencies = T)
library(gage)
library(pathview)


## Find species in KEGG database
# org <- "Homo sapiens"
# species <- unlist(sapply(1:ncol(korg), function(i) {
#   agrep(org, korg[, i])
# }))
# korg[species, 1, drop = F]
# 
# kegg_arab <- kegg.gsets("hsa")
# kegg_human <- kegg_arab$kg.sets[kegg_arab$sigmet.idx]
# save(kegg_human, file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/kegg_human.RData")
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/kegg_human.rdata")
load("~/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/x.rdata")

## Make count database
rnaseq_counts <- as.data.frame(x$counts)

# Remove Zero counts
non_zero <- rowSums(rnaseq_counts) != 0
rnaseq_counts <- rnaseq_counts[non_zero,]

# Calculate library sizes and normalise
libsizes <- colSums(rnaseq_counts)
size_factor <- libsizes / exp(mean(log(libsizes)))
norm_counts <- t(t(rnaseq_counts) / size_factor)

# Stabilise variance (add 8 as per GAGE manual to avoid -inf)
norm_counts <- log2(norm_counts + 8)

## Pick the sampes you want to differentially compare
ref_idx <- grep("CD27HI", colnames(norm_counts)) #CD27HI (naive)
samp_idx <- grep("CD27LO", colnames(norm_counts)) #CD27LO (effectors) 

# Enrichment analysis
native_kegg_fc <- gage(norm_counts, 
                       gsets = kegg_human, 
                       ref = ref_idx, 
                       samp = samp_idx, 
                       compare ="unpaired")

## First four pathways which are downregulated (up in Naive)
head(native_kegg_fc$less[,1:5], 4)

# Significance test to find pathways that are up and down-regulated
wt_hy_sig_kegg <- sigGeneSet(native_kegg_fc, outname = "wt_hy.kegg")


log_fc = norm_counts[, samp_idx]-rowMeans(norm_counts[, ref_idx])


greater_set <- native_kegg_fc$greater[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$greater[, "q.val"])
greater_ids <- rownames(native_kegg_fc$greater)[greater_set]
head(greater_ids)

less_set <- native_kegg_fc$less[, "q.val"] < 0.1 &
  !is.na(native_kegg_fc$less[,"q.val"])
less_ids <- rownames(native_kegg_fc$less)[less_set]
head(less_ids)

combine_ids <- substr(c(greater_ids, less_ids), 1, 8)
greater_ids1 <- substr(greater_ids, 1, 8)
less_ids1 <- substr(less_ids, 1, 8)


browseVignettes("gage")
# Bringing in EdgeR
d <- x

# Estimate the dispersion
d1 <- estimateCommonDisp(d, verbose = T)
d1 <- estimateTagwiseDisp(d1)
# plotBCV(d1)

# Differential gene expression
de.com <- exactTest(d1, pair = c("VD1.CD27HI", "VD1.CD27LO"))
topTags(de.com, n = 5)

FDR <- p.adjust(de.com$table$PValue, method = "BH")
de1 <- decideTestsDGE(de.com, adjust.method = "BH", p.value = 0.05)
summary(de1)

# Make suitable for pathway analysis
isDE <- as.logical(de1) # covert to DE set to true/false set

DEnames <- rownames(d)[isDE] # get the DE gene names
edger.fc <- de.com$table$logFC # get the log fold change
names(edger.fc) <- rownames(de.com$table) # assign row names to fold change 
exp.fc <- edger.fc
head(exp.fc)

DE_foldchange <- edger.fc[DEnames]

length(DE_foldchange)

# GAGE Analysis
fc.kegg.p <- gage(exp.fc, gsets = kegg_human, ref = NULL, samp = NULL)
head(fc.kegg.p$greater[,1:5], 4)
wt_hy_sig_kegg <-sigGeneSet(fc.kegg.p, outname = "wt_hy.kegg")


setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Figures/Geneset_Enrichment")
# pv.out.list <- sapply(combine_ids[1:3], function(pid) pathview(
#   gene.data =  exp.fc, pathway.id = pid,
#   species = "hsa", out.suffix = "27LO_compared_to_27HI",
#   gene.idtype = "KEGG")) # Gives loads more shit too...

pv_replication <- pathview(gene.data = DE_foldchange, 
                           gene.idtype = "KEGG", 
                           pathway.id = combine_ids[1:3], 
                           species = "hsa", 
                           out.suffix = "CD27LO_compred_to_CD27HI", 
                           keys.align = "y", 
                           kegg.native = T, 
                           match.data = T, 
                           key.pos = "topright")

