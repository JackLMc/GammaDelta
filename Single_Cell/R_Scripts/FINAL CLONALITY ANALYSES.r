
library(Seurat)

load("Final_seurat.RData")

clonality_values = read.table("../metadata//cell_clonality_from_iRep.csv", sep=",", row.names = 1, header=TRUE)

head(clonality_values)

data_noBatch <- AddMetaData(data_noBatch, clonality_values[rownames(data_noBatch@data.info),])

head(data_noBatch@data.info)

FeaturePlot(data_noBatch, c("D_clonality", 'G_clonality', 'max_clonality', 'mean_clonality'), reduction.use='pca')

EMRA_clonality = data_noBatch@data.info[data_noBatch@data.info$sort_class =='EMRA',]

head(EMRA_clonality)

hist(EMRA_clonality$D_clonality)

hist(EMRA_clonality$G_clonality)

hist(EMRA_clonality$max_clonality)

hist(EMRA_clonality$mean_clonality)

classify_clonality <- function(x, low_ceiling=10, high_floor=20){
    max_clonality = as.numeric(x['max_clonality'])
    sort_class = x['sort_class']
    clone_class = 'no_class'
    #if (max_clonality>0 && sort_class == 'EMRA') {
    if (max_clonality>0) {
        if (max_clonality < low_ceiling){
            clone_class = 'low_clonality'
        } else if (max_clonality>high_floor){
            clone_class = 'high_clonality'
        } else{
        } 
    }
    return(clone_class)
}

data_noBatch <- AddMetaData(data_noBatch, apply(data_noBatch@data.info, 1, classify_clonality, low_ceiling=20, high_floor=20), "clone_class_all")

head(data_noBatch@data.info)

classify_clonality <- function(x, low_ceiling=10, high_floor=20){
    max_clonality = as.numeric(x['max_clonality'])
    sort_class = x['sort_class']
    clone_class = 'no_class'
    if (max_clonality>0 && sort_class == 'EMRA') {
        if (max_clonality < low_ceiling){
            clone_class = 'low_clonality'
        } else if (max_clonality>high_floor){
            clone_class = 'high_clonality'
        } else{
        } 
    }
    return(clone_class)
}

data_noBatch <- AddMetaData(data_noBatch, apply(data_noBatch@data.info, 1, classify_clonality, low_ceiling=20, high_floor=20), "clone_class_EMRA")

PCAPlot(data_noBatch, do.return=T, group.by="clone_class_all", no.legend=F)

PCAPlot(data_noBatch, do.return=T, group.by="clone_class_EMRA", no.legend=F)

data_noBatch <- AddMetaData(data_noBatch, apply(data_noBatch@data.info, 1, classify_clonality, low_ceiling=10, high_floor=20), "clone_class_EMRA")

PCAPlot(data_noBatch, do.return=T, group.by="clone_class_EMRA", no.legend=F)

data_noBatch = SetAllIdent(data_noBatch, "clone_class_EMRA")

high_clonality_markers = FindMarkers(data_noBatch, 'high_clonality', 'low_clonality', min.pct=0.25, only.pos=T)

low_clonality_markers = FindMarkers(data_noBatch, 'low_clonality', 'high_clonality', min.pct=0.25, only.pos=T)

write.table(high_clonality_markers, file='high_clonality_markers.txt', sep="\t", quote=F)

write.table(low_clonality_markers, file='low_clonality_markers.txt', sep="\t", quote=F)

VlnPlot(data_noBatch, c('PRF1', 'GZMA', 'CX3CR1', 'SELL', 'TCF7', 'LEF1', 'IL7R', 'LTB', 'GZMB'), nCol=2)

high_clonality_markers[c('PRF1', 'GZMA', 'GZMB', 'CX3CR1', 'TCF7'), ]

clonotype_metadata = read.table("../metadata/clonotype_metadata.txt", row.names = 1)

colnames(clonotype_metadata) = 'clonotype'

data_noBatch <- AddMetaData(data_noBatch, clonotype_metadata[rownames(data_noBatch@data.info),,drop=FALSE], 'clonotype')

head(data_noBatch@data.info)

test = PCAPlot(data_noBatch, group.by="clonotype", do.return=T)
test2 = ggplot_build(test)

fill_colours = names(sort(table(test2$data[[1]]['colour']), decreasing=T))[2:14]

cols = append(fill_colours, c("#D3D3D3"))

PCAPlot(data_noBatch, group.by="clonotype", do.return=T, cols.use=cols)

table(data_noBatch@data.info['clonotype'])

data_noBatch = SetAllIdent(data_noBatch, "clonotype")

clonotype_A_markers = FindMarkers(data_noBatch, 'A', c("E","J","K"))

clonotype_E_markers = FindMarkers(data_noBatch, 'E', c("A","J","K"))

clonotype_J_markers = FindMarkers(data_noBatch, 'J', c("A","E","K"))

clonotype_K_markers = FindMarkers(data_noBatch, 'K', c("A","E","K"))

clonotype_AE_markers = FindMarkers(data_noBatch, c('A','E'), c("J","K"))

clonotype_AJ_markers = FindMarkers(data_noBatch, c('A','J'), c("E","K"))

clonotype_AK_markers = FindMarkers(data_noBatch, c('A','K'), c("E","J"))

write.table(clonotype_A_markers, 'clonotype_A_markers.txt', sep="\t", quote=F)
write.table(clonotype_E_markers, 'clonotype_E_markers.txt', sep="\t", quote=F)
write.table(clonotype_J_markers, 'clonotype_J_markers.txt', sep="\t", quote=F)
write.table(clonotype_K_markers, 'clonotype_K_markers.txt', sep="\t", quote=F)
write.table(clonotype_AE_markers, 'clonotype_AE_markers.txt', sep="\t", quote=F)
write.table(clonotype_AJ_markers, 'clonotype_AJmarkers.txt', sep="\t", quote=F)
write.table(clonotype_AK_markers, 'clonotype_AK_markers.txt', sep="\t", quote=F)

genes_for_violin = list(NK_receptors=c("KLRB1","KLRG1","KLRC2","KLRC3","KLRC4","CADM1","FCGR3A","FCGR3B","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DS2","KIR3DS1","KIR3DL2","CD160","KLRF1","KLRK1","KLRD1","KLRC1","LILRB1","LILRB2"), cytokines_chemokines_and_receptors=c("CCL28","CCL3","CCL4","CCL5","TNFSF13B","TGFBR1","IL32","CXCR4","IL7R","CX3CR1"), homing_adhesion=c("ITGA5","ITGAD","ITGAX","ITGAE","ITGAL","SELL","PECAM1","S1PR5","CD2","SEMA4C","SEMA4F","SEMA3C","CD52","CLEC11A"), effector_molecules=c("PRF1","GZMA","GZMB","GZMK","GZMM","GZMH"), transcription=c("NANOG","TCF7","PRDM1","EOMES","SOCS2"))

genes_for_violin

pdf(file="test.pdf")
VlnPlot(data_noBatch,c("KLRB1", "KLRG1", "KLRC2", "KLRC3"), size.x.use=8, size.y.use=8, size.title.use=8, size.use=0.75, y.log=F, group.by = "clonotype")
VlnPlot(data_noBatch,c("KLRC4", "CADM1", "FCGR3A", "FCGR3B"), size.x.use=8, size.y.use=8, size.title.use=8, size.use=0.75, y.log=F, group.by = "clonotype")
dev.off()

test = genes_for_violin[['NK_receptors']]

test

test = split(test, ceiling(seq_along(test)/4))

vp = function(g, data_to_use=data_noBatch){
    VlnPlot(data_to_use, g, size.x.use=8, size.y.use=8, size.title.use=8, size.use=0.75, y.log=F, group.by = "clonotype")
}

vp(test[[1]])

save.image("final_clonality_analysis.RData")

load("final_clonality_analysis.RData")

for (n in names(genes_for_violin)){ 
    file = paste("violin_plots/", n, ".pdf", sep="")
    pdf(file=file)
    genes = genes_for_violin[[n]]
    genes = split(genes, ceiling(seq_along(genes)/4))
    lapply(genes, vp)
    dev.off()
}


VlnPlot(data_noBatch, c('BACH2', 'IFNG', 'ID2', 'ID3', 'BATF', 'IRF4'),  nCol=2)

FeaturePlot(data_noBatch, c('ID2'), reduction.use='pca')

tpms['IFNG',]
