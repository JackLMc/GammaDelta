library(tidyverse)
library(ggpubr)
load("Bulk/Final_edgeR.RData")

# Violins of important stuff figure 3A
counts_ <- as.data.frame(x$counts) %>% rownames_to_column(., var = "entrezgene_id")
# 
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- counts_$entrezgene_id

# listAttributes(mart)[grep("entrez", listAttributes(mart)$name), ]

G_list <- getBM(filters = "entrezgene_id", attributes = c("entrezgene_id", "external_gene_name"),
                values = genes, mart = mart)

# backup_g_list <- G_list
# G_list <- backup_g_list
# counts1 <- merge(counts_, G_list, by = "entrezgene_id")
# 
# head(counts1)
# 
# rescale <- function(x) {
#   s = sum(x)
#   (x/s) * 1e6
# }
# 
# 
# head(G_list)
# head(counts_)
# 
# counts_ <- apply(counts_, 2, rescale)
# 

cpms <- as.data.frame(cpm) %>% rownames_to_column(., var = "entrezgene_id")

cpms1 <- merge(cpms, G_list, by = "entrezgene_id") %>% 
  gather(key = "Sample", value = "cpm", -entrezgene_id, -external_gene_name)
cpms1$log_cpm <- log(cpms1$cpm + 1)

samp_data <- as.data.frame(x$samples) %>% rownames_to_column(., "Sample")
samp_data <- samp_data[, c("Sample", "group")]

cpms2 <- merge(cpms1, samp_data, by = "Sample")

# cpms3[grep("CCL5", cpms3$external_gene_name),]

cpms3 <- droplevels(subset(cpms2, group != "VD2"))
cpms3$Cell <- ifelse((grepl("CD8", cpms3$Sample)), "CD8", "VD1")

dodge = position_dodge(width = 0.9)

cpms3$group <- factor(cpms3$group,levels = c("CD8.Naive", "CD8.EMRA", "VD1.CD27HI", "VD1.CD27LO"))


these_genes <- c("TCF7", "LEF1", "MYC", "TBX21", "EOMES", "PRDM1")

# cpms3[grep("TNF", cpms3$external_gene_name),]

i <- 1
for(i in 1:length(these_genes)){
  j <- these_genes[i]
  working_df <- droplevels(subset(cpms3, external_gene_name == j))
  dodge <- position_dodge(width = 0.9)
  temp_plot <- ggplot(working_df, aes(x = Cell, y = log_cpm, fill = group)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.5) +
    stat_summary(fun.y = mean, geom = "bar", position = dodge) + 
    theme_bw() + 
    scale_fill_manual(values = cbcols) +
    # labs(x = "Tissue", y = "FPKM") + 
    stat_compare_means(label = "p.signif", label.y = max(working_df$log_cpm + 0.5)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "none", text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(working_df$log_cpm + 1))) + 
    labs(x = "Cell population", y = paste0("log(cpm + 1) of ", j))
  filen <- paste0(j, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Bulk/Figures/CPMs/",
                   height = 6, width = 6, units = "in")
}

# For the legend
pdf("../Figures/CPMs/Legend.pdf")
Cd27_df <- droplevels(subset(cpms3, external_gene_name == "CD27"))

ggplot(Cd27_df, aes(x = Cell, y = log_cpm, fill = group)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.5) +
  stat_summary(fun.y = mean, geom = "bar", position = dodge) + 
  theme_bw() + 
  scale_fill_manual(values = cbcols) +
  # labs(x = "Tissue", y = "FPKM") + 
  stat_compare_means(label = "p.signif", label.y = max(Cd27_df$log_cpm + 0.5), method = "wilcox.test") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "top", text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(Cd27_df$log_cpm + 1))) + 
  labs(x = "Cell population", y = paste0("log(cpm + 1) of ", "CD27"))



ggplot(droplevels(subset(cpms3, external_gene_name == "CD27")), aes(x = Cell, y = log_cpm, fill = group)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.5) +
  stat_summary(fun.y = mean, geom = "bar", position = dodge) + 
  scale_fill_manual(values = cbcols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top", text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 9)) 
dev.off()



# All the cytos
these_cytos <- c("IFNG", "TNF", "IL10", "CCL5")

# cpms3[grep("TNF", cpms3$external_gene_name),]


i <- 1
for(i in 1:length(these_cytos)){
  j <- these_cytos[i]
  working_df <- droplevels(subset(cpms3, external_gene_name == j))
  dodge <- position_dodge(width = 0.9)
  temp_plot <- ggplot(working_df, aes(x = Cell, y = log_cpm, fill = group)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.5) +
    stat_summary(fun.y = mean, geom = "bar", position = dodge) + 
    theme_bw() + 
    scale_fill_manual(values = cbcols) +
    # labs(x = "Tissue", y = "FPKM") + 
    stat_compare_means(label = "p.signif", label.y = max(working_df$log_cpm + 0.5)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "none", text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(working_df$log_cpm + 1))) + 
    labs(x = "Cell population", y = paste0("log(cpm + 1) of ", j))
  filen <- paste0(j, ".pdf")
    ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                     path = "./Bulk/Figures/CPMs/Cytokines/",
                     height = 6, width = 6, units = "in")
  }



