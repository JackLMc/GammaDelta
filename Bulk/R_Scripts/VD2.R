load("./Bulk/Final_edgeR.RData")

head(CD27LO.vs.VD2)

inte <- droplevels(subset(CD27HI.vs.Naive, adj.P.Val <= 0.05))
View(inte)
for_plot <- droplevels(subset(CD27LO.vs.VD2, adj.P.Val <= 0.05))
this <- droplevels(subset(CD27HI.vs.VD2, adj.P.Val <= 0.05))
View(this)


for_plo <- droplevels(subset(EMRA.vs.VD2, adj.P.Val <= 0.05))

head(for_plot)
View(for_plo)

for_plot$logFC <- for_plot$logFC * -1
head(for_plot)


ord_lev <- for_plot$SYMBOL[order(for_plot$adj.P.Val)]

for_plot$SYMBOL <- factor(for_plot$SYMBOL, levels = ord_lev)

library(cowplot)

pdf("../Figures/VD2/Pval_VD2_vs_CD27LO.pdf", height = 10, width = 14)
ggplot(for_plot, aes(x = SYMBOL, y = logFC)) +
  geom_bar(stat = "identity", width = 0.95) + 
  # scale_fill_manual(values = cbcols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", text = element_text(size = 10)) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  panel_border(remove = T) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 0) +
  geom_label(label = levels(for_plot$SYMBOL), nudge_y = 5) + 
  ylim(-5, 17)
dev.off()




## subset to only contain genes in Common_genes
de.common <- which(dt[, "CD27LOvsCD27HI"]!=0 & dt[,"EMRAvsNaive"]!=0)
Common_genes <- tfit$genes$ENTREZID[de.common]
i <- which(v$genes$ENTREZID %in% Common_genes)

## Define needed bits
library(gplots)
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell1 <- c("#56B4E9", "#E69F00", "#009E73", "#999999", "#CC79A7")[v$targets$group]
data.frame(v$targets$group, col.cell1)
mycol <- colorpanel(1000,"blue","white","red")

pdf("../Figures/VD2/426_genes_VD2.pdf", width = 10, height = 50)
heatmap.2(v$E[i,], scale = "row",
          labRow = ifelse((is.na(v$genes$SYMBOL[i])), v$genes$ENTREZID[i], v$genes$SYMBOL[i]), labCol = NA,
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8, 6), lhei = c(2, 10), 
          hclustfun = hclustAvg, ColSideColors = col.cell1)
dev.off()

##### VD2 Stuff
## Hierarchical clustering - no longer produces nice one
library(tidyverse)
library(ggdendro)

v$E

head(DF)
DF <- as.data.frame(v$E)
length(Common_genes)

DF1 <- DF[row.names(DF) %in% Common_genes, ] %>%
  rownames_to_column(var = "Genes") %>% 
  gather(contains("."), key = "ID" , value = "value") %>%
  spread(key = "Genes", value = "value") 


head(DF1)
DF2 <- data.frame(DF1[, names(DF1) != "ID"], row.names = DF1[, names(DF1) == "ID"])
hc <- hclust(dist(DF2), "average")
?hclust
p1 <- ggdendrogram(hc, rotate = F, size = 2, leaf_labels = F)

df2 <- data.frame(cluster = cutree(hc, 1), cell.types = factor(hc$labels, levels = hc$labels[hc$order]))

merging <- as.data.frame(x$samples) 
merging1 <- rownames_to_column(merging, var = "cell.types")

merging2 <- merging1[, c("cell.types", "group")]

df3 <- merge(df2, merging2, by = "cell.types")

p2 <- ggplot(df3, aes(cell.types, y = 1, fill = factor(cluster))) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), 
                    name = "Cluster",
                    breaks = c("1", "2"),
                    labels = c("Cluster 1", "Cluster 2"))
p3 <- ggplot(df3, aes(cell.types, y = 1, fill = factor(group))) + 
  geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = cbcols, 
                    name = "Cluster")
dev.off()
p3
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)  
gp3 <- ggplotGrob(p3)

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

library(gridExtra)
pdf("../../Bulk/Figures/VD2/Dendro.pdf")
grid.arrange(gp1, gp3, ncol = 1, heights = c(4/5, 1/5, 1/5)) # Doesn't 
dev.off()






