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


load("../Final_edgeR.RData")

## Showing the removal of the lowly expressed transcripts
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow = c(1,2))
plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2,
     main = "", xlab = "")
title(main = "A. Raw data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", samplenames, text.col = col, bty = "n")
lcpm <- cpm(x, log = T)
plot(density(lcpm[,1]), col = col[1], lwd = 2, ylim = c(0,0.4), las = 2,
     main = "", xlab = "")
title(main = "B. Filtered data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", samplenames, text.col = col, bty = "n")


## Showing theoretical effect of normalisation
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05) # reduce first sample to 5%
x2$counts[,2] <- x2$counts[,2]*5 # Inflate second sample by x5

### Graph
par(mfrow = c(1,2))
lcpm <- cpm(x2, log = T)
boxplot(lcpm, las = 2, col = col, main = "")
title(main = "A. Example: Unnormalised data", ylab = "Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors
lcpm <- cpm(x2, log = T)
boxplot(lcpm, las = 2, col = col, main = "")
title(main = "B. Example: Normalised data", ylab = "Log-cpm")


# Unsupervised clustering of samples & runs
library(RColorBrewer)
par(mfrow = c(1,2))
col.group <- as.factor(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- as.factor(lane)
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels = group, col = col.group)
title(main = "A. Sample groups")
plotMDS(lcpm, labels = lane, col = col.lane, dim = c(1,2))
title(main = "B. Sequencing lanes")

## Online launch of this.
# biocLite("Glimma", dependencies = T)
library(Glimma)
# glMDSPlot(lcpm, labels = paste(group, lane, sep = "_"),
#           groups = x$samples[,c(2,5)], launch = F)

col.cell <- c("#999999","#56B4E9","#E69F00","#009E73","#CC79A7")[x$samples$group]
data.frame(x$samples$group, col.cell)

png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Bulk/Figures/Paper/PCA_of_all_genes.png", width = 150, height = 150, units = "mm", res = 300)
plotMDS(lcpm, pch = 16, cex = 1.2, col = col.cell, top = 11999)
legend("top",
       fill = c("#999999", "#56B4E9",
                "#E69F00", "#009E73",
                "#CC79A7"),
       legend = levels(x$samples$group))
# Add a title
title(main = "A. MDS of log(cpm) for all genes")
dev.off()
