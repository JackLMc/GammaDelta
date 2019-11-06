## INSTRUCTIONS
# Export a file from flowjo (scale values I believe) with the name in the format Biologicalcondition_timepoint_individual.csv
path_to_files <- "./Martin/Data/CSV/Panel1/Data"

# Set up a directories
getwd() # make sure this matches where you think you are i.e. the place where you're exporting to.

# RUN ONCE
# setwd(paste0(path_to_files, "/"))
# 
# subDir <- c("Data", "Figures")
# for(i in 1:length(subDir)){
#   j <- subDir[i]
#   if (!file.exists(j)){
#     dir.create(file.path(getwd(), j))
#   }}
# 
# setwd("../../../../")

# Load libraries as you go
library(devtools)
# install_github("tchitchek-lab/SPADEVizR")
library(SPADEVizR)
library(tidyverse)

# Read in samples into a list, take name of the file as the list name
filelist <- list.files(path = path_to_files, pattern = "csv$", full.names = T) %>% unlist

# Read in files and combine
lists <- lapply(filelist, read.csv, header = T, stringsAsFactors = F)
names(lists) <- filelist
names(lists) <- gsub(paste0(path_to_files, "/"), "", names(lists)) 
names(lists) <- gsub(".csv", "", names(lists))

lists1 <- lapply(names(lists), function(n, x){
  x[[n]]$Sample <- n
  return (x[[n]])},
  lists)

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)}

full_data <- multi_join(lists1, full_join)
remove_cols <- c("Time", "FSC.A", "FSC.H", "SSC.A", "SSC.H", "Comp.BV510.A....DEAD", "Comp.PE.Cy7.A....VD1") # VD1 not actually included
subs_data <- full_data[, !('%in%'(colnames(full_data), remove_cols))]

colours <- colnames(subs_data)
colnames(subs_data) <- gsub(".*\\.", "", colnames(subs_data))

# Arcsin transformations
library(flowCore)
ncol(subs_data)

asinhframe <- list()

colours

asinhframe[["VD1"]]  <- c(a = 1, b = 0.005, c = 0) #FL1 (FITC/Alexa.Fluor.488)
asinhframe[["VD2"]]  <- c(a = 1, b = 0.005, c = 0) #FL2 (PE)
asinhframe[["CD27"]]  <- c(a = 1, b = 0.002, c = 0) #FL3 (Texas Red?)
# asinhframe[[8]]  <- c(a = 1, b = 0.002, c = 0) #FL4 (PerCP-Cy5.5)
# asinhframe[["VD1"]]  <- c(a = 1, b = 0.002, c = 0) #FL5 (PE Cy7)
asinhframe[["eomes"]] <- c(a = 1, b = 0.002, c = 0) # FL6 (APC/AF647)
# asinhframe[[11]] <- c(a = 1, b = 0.001, c = 0) # FL7 (AF680/AF700)
asinhframe[["TIGIT"]] <- c(a = 1, b = 0.002, c = 0) # FL3 (APC Cy 7)
asinhframe[["CD45RO"]] <- c(a = 1, b = 0.002, c = 0) #FL9 (BV421)
# asinhframe[[14]] <- c(a = 1, b = 0.001, c = 0) #FL10 (BV510)


asinhframe[["IL7Ra"]] <- c(a = 1, b = 0.001, c = 0) # (BV605)
asinhframe[["CD8"]] <- c(a = 1, b = 0.001, c = 0) # (BV650)
asinhframe[["CD45RA"]] <- c(a = 1, b = 0.001, c = 0) # (BV711)
asinhframe[["CD3"]] <- c(a = 1, b = 0.001, c = 0) # (BV786)

asinhframe[["Vg9"]] <- c(a = 1, b = 0.001, c = 0) # (PE.Cy5)

missing <- c("Comp.BV605.A....IL7Ra", "Comp.BV650.A....CD8", "Comp.BV711.A....CD45RA",
             "Comp.BV786.A....CD3", "Comp.PE.Cy5.A....Vg9") # guessed


df <- subs_data

for (x in colnames(df)[unlist(lapply(df, is.numeric))]){
  aSinh <- arcsinhTransform(a = asinhframe[[x]][1],
                            b = asinhframe[[x]][2],
                            c = asinhframe[[x]][3])
  df[,x] <- aSinh(df[,x])
  print(x)}

library(ggpubr)
ggdensity(df$eomes)
ggdensity(df$TIGIT)
ggdensity(df$VD1)
ggdensity(df$CD45RO)
ggdensity(df$IL7Ra)
ggdensity(df$CD8)
ggdensity(df$CD45RA)
ggdensity(df$CD3)
ggdensity(df$VD2)
ggdensity(df$CD27)
ggdensity(df$Vg9)

range(df$CD27) # Needs tightening


plot(df[,"VD1"],
     df[, "CD27"])

plot(df[,"eomes"],
     df[,"TIGIT"])

plot(df[,"VD1"],
     df[,"VD2"])

plot(df[,"Vg9"],
     df[,"VD2"])

plot(df[,"CD3"],
     df[,"VD2"])



# Create a unique ID 
## colnames(subs_data)[colnames(subs_data) == "SampleID"] <- "CellID"
## subs_data$CellID <- paste("X", subs_data$CellID, sep = "")
## subs_data1 <- column_to_rownames(subs_data, var = "CellID") ## for some reason the SampleID that comes with is not unique? Unsure why - just create a unique one
data_clean <- df

rownames(data_clean) <- paste0("X", rownames(data_clean))

Cells_and_samples <- rownames_to_column(data_clean, var = "CellID")
CaS <- Cells_and_samples[, c("CellID", "Sample")]

## NEED TO BOOLEAN GATE, BASICALLY THE CELLS WHICH ARE IN BOTH YOU NEED TO DROP THE ONES WHICH ARE IN THE CD8... THESE ARE ACTIVATED GD
# DUPLICATED DOES THIS, AS IT GETS RID OF THE SECOND INSTANCE OF IT BUT THIS IS UNSAFE

# Clustering
library(Rphenograph)
set.seed(1)

a <- cbind(data_clean, cluster = factor(Rphenograph(as.matrix(data_clean[, !('%in%'(colnames(data_clean), "Sample"))]))[[2]]$membership))
a$cluster <- paste0("cluster", a$cluster) %>% as.factor()

# RTsne
head(data_clean)
# data_clean[duplicated(data_clean[, !('%in%'(colnames(data_clean), "Sample"))]),]
# duplicated(data_clean) | duplicated(data_clean[nrow(df):1, ])[nrow(data_clean):1]


library(Rtsne)
tsne_out <- Rtsne(as.matrix(as.matrix(data_clean[, !('%in%'(colnames(data_clean), "Sample"))]))) # Don't have to run this... v. slow
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")
Pcluster <- a[, "cluster"]

# tSNE plot
pdf("./Desktop/FCS/Figures/tSNE.pdf") # CHANGE THIS TO A PLACE YOUR WANT TO SAVE YOUR FIGURES
ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()



marker <- a[, "CD27"]

ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = marker)) + 
  geom_point() + 
  scale_colour_gradient(low = "#E3E3E3", high = "#413CFF",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(axis.line.x = element_line(color = "black", size = .3),
        axis.line.y = element_line(color = "black", size = .3))








# Making the right dataframe for import to SPADEVizR
## We'll be using the importResultsFromTables() function
data_clean1 <- a
data_clean1$Sample_cluster <- paste(data_clean1$Sample, data_clean1$cluster, sep = "___") %>% as.factor()
head(data_clean1$Sample_cluster)

Cell_numbers <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(data_clean1$Sample_cluster)){
  print(paste0("Working on: ", i))
  work <- droplevels(subset(data_clean1, Sample_cluster == i))
  number_of_cells <- nrow(work)
  Cell_numbers[c, "ID"] <- i
  Cell_numbers[c, "Number"] <- number_of_cells
  c <- c + 1
}

cluster.abundances <- separate(Cell_numbers, "ID", into = c("Sample", "cluster"), sep = "___") %>% 
  spread(key = "Sample", value = "Number") %>%
  column_to_rownames(var = "cluster")


# all.equal(unname(sapply(data_clean1[, unlist(lapply(data_clean1, is.numeric))], median, na.rm = T)["Yb171Di....171Yb_DNAM"]), 
#           median(data_clean$Yb171Di....171Yb_DNAM)) ## Check
# columns <- paste(paste0(colnames(data_clean1[, unlist(lapply(data_clean1, is.numeric))]), " = double()"), collapse = ", ") # Doesn't work too nicely

pheno_list <- list()
c <- 1
for(i in levels(data_clean1$Sample_cluster)){
  work <- droplevels(subset(data_clean1, Sample_cluster == i))
  medians_ <- sapply(work[, unlist(lapply(work, is.numeric))], median, na.rm = T)
  pheno_list[[i]] <- as.data.frame(medians_)
  c <- c + 1
}


pheno_list1 <- lapply(seq_along(pheno_list), function(i) {
  colnames(pheno_list[[i]]) <- names(pheno_list)[i]
  pheno_list[[i]]})


cluster.phenotypes <- do.call("cbind", pheno_list1) %>% t() %>% as.data.frame() %>% rownames_to_column(., var = "ID") %>%
  separate(., col = "ID", into = c("sample", "cluster"), sep = "___")

write.table(cluster.abundances, file = "./Desktop/FCS/Data/cluster.abundances.txt", sep = "\t") # CHANGE THIS TO A PLACE YOU WANT TO SAVE YOUR DATA FILES.
write.table(cluster.phenotypes, file = "./Desktop/FCS/Data/cluster.phenotypes.txt", sep = "\t")


# START THE SPADE BIT
library(SPADEVizR)
results_tables <- importResultsFromTables(cluster.abundances = cluster.abundances, cluster.phenotypes = cluster.phenotypes)

## Add some stuff to do with them (name file by biological condition, timepoint, individuals, just important stuff to do with that)
assignments <- data.frame(ID = c(colnames(cluster.abundances)))
assignments <- separate(assignments, "ID", into = c("bc", "tp", "ind"), sep = "_") # Very annoying they only take these columns.. But hey ho
rownames(assignments) <- colnames(cluster.abundances)

results <- assignContext(results_tables, assignments = assignments)

# I CAN'T DO ANYMORE, CAUSE I ONLY HAVE ONE SAMPLE... SHOULD WORK THOUGH. JUST FOLLOW THE ONLINE STUFF...
samples   <- c(colnames(cluster.abundances))
resultsAC <- identifyAC(results, samples = samples, mu = 1, th.pvalue = 0.01) 

samples <- c(colnames(cluster.abundances))
treeViewer(results, samples = samples, highlight = resultsDAC, marker = "Yb172Di....172Yb_CX3CR1")

