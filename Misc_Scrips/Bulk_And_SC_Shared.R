getwd()

SC <- read.csv("./Single_Cell/Output/SC_DE.csv")
Bulk <- read.csv("./Bulk/Output/CD27LO.vs.CD27HI.csv")

Bulk <- droplevels(subset(Bulk, adj.P.Val <= 0.05))
Bulk_and_SC_Shared_genes <- Bulk[Bulk$SYMBOL %in% as.character(droplevels(Bulk$SYMBOL[Bulk$SYMBOL %in% SC$SYMBOL])),]



View(Bulk_and_SC_Shared_genes)
write.csv("Bulk_and_SC_Shared_genes.csv", x = as.data.frame(Bulk_and_SC_Shared_genes), row.names = F)

