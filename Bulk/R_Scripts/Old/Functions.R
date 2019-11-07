# FUNCTIONS Bioportal
# Install and load required packages
required <- c("tidyverse",
              #"fields",
              "reshape2",
              "ggpubr",
              "data.table",
              "gplots",
              "devtools",
              "randomForest",
              "factoextra",
              "Rsubread",
              "ggbiplot",
              "Rsubread")
for (lib in required)
{
  if (!require(lib, character.only = T))
  {
    install.packages(lib, dependencies = T)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

# Comparisons
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"),
                       c("VD1.CD27HI", "CD8.EMRA"),
                       c("VD1.CD27HI", "CD8.Naive"),
                       c("VD1.CD27HI", "VD2"),
                       c("VD1.CD27LO", "CD8.EMRA"),
                       c("VD1.CD27LO", "CD8.Naive"),
                       c("VD1.CD27LO", "VD2"),
                       c("CD8.EMRA", "CD8.Naive"),
                       c("CD8.EMRA", "VD2"),
                       c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999",
            "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00",
            "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")


produce_pca <- function(num){
  df <- as.data.frame(num)
  df1 <- rownames_to_column(df, var = "gene")
  df2 <- data.frame(t(df1[-1]))
  colnames(df2) <- df1[, 1]
  df3 <- rownames_to_column(df2, var = "Temp")
  df3$group <- ifelse(grepl("VD1.CD27LO", df3$Temp) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD27LO", df3$Temp), "VD1.CD27LO",
                  ifelse(grepl("CD8.EMRA", df3$Temp) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{4}[[:punct:]]{1}EMRA", df3$Temp), "CD8.EMRA",
                         ifelse(grepl("CD8.Naive", df3$Temp) | grepl("CD8[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]Naive", df3$Temp), "CD8.Naive",
                                ifelse(grepl("VD1.CD27HI", df3$Temp) | grepl("VD1[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]{1}[[:digit:]]{2}[[:punct:]]CD2HI", df3$Temp), "VD1.CD27HI",
                                       ifelse(grepl("VD2", df3$Temp), "VD2", "none")))))
  df3a <- data.frame(df3[, names(df3) != "Temp"], row.names = df3[, names(df3) == "Temp"])
  prin_comp <- prcomp(df3a[, names(df3a) != "group"], scale. = T)
  sgroup <- df3a[, "group"]
  g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
                groups = sgroup, ellipse = T,
                circle = T,
                var.axes = F
  )
  g <- g + scale_color_manual(values = cbcols)
  g <- g + theme_bw()
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  
  g <- g + ggtitle("PCA")
  return(g)
}
# cbPalette <- c("#999999" grey, "#E69F00" orange, "#56B4E9" light blue, "#009E73" dark green, "#F0E442" yellow,
# "#0072B2" navy, "#D55E00" red, "#CC79A7" magenta)


# Get the Bam files
get_sorbam <- function(foldertobam){
  thousand.folders <- list.dirs(path = foldertobam, full.names = T)
  filelist1 <- sapply(thousand.folders[-1], function(x)
  {list.files(x, pattern = ".sor.bam$", full.names = T)
  })
  bam.files <- unlist(filelist1)
  setwd("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/GammaDelta/Bulk")
  return(bam.files)
}

# Take only significant values at the correct fold change
Take_Sigs <- function(df){
  try <- droplevels(subset(df, adj.P.Val < 0.05 & logFC > log2(2.5) |
                             adj.P.Val < 0.05 & logFC < log2(2.5) | 
                             adj.P.Val < 0.05 & logFC == log2(2.5) |
                             adj.P.Val < 0.05 & logFC == log2(2.5) |
                             adj.P.Val == 0.05 & logFC > log2(2.5) |
                             adj.P.Val == 0.05 & logFC < log2(2.5) | 
                             adj.P.Val == 0.05 & logFC == log2(2.5) |
                             adj.P.Val == 0.05 & logFC == log2(2.5)))
  return(try)}

# Trim whitespace
trimWS <- function (x) gsub("^\\s+|\\s+$", "", x)

# Write a CSV for Raw data (with space titles for raw data)
writeCsvD <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("C:/Users/User/OneDrive/University of Birmingham/PhD/TCGA/Data/",fn,".csv"), row.names = F)}

# Write a csv for non raw data
writeCsvO <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/GammaDelta/Bulk/Output/DEgenes/",fn,".csv"), row.names = F)}

writeCsvA <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/Desktop/Attempt/grch37/Sorted/",fn,".csv"), row.names = F)}

writeCsvCount <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/GammaDelta/Bulk/Counts/",fn,".csv"), row.names = F)}

writeCsvGS <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/GammaDelta/Bulk/Output/GSEA/",fn,".csv"), row.names = F)}

# Rename column names
renamecolumn <- function(df, oldname, newname){
  names(df)[names(df) == oldname] <- newname
  return(df)
}

# Function to remove the addition of X. column at the start of a dataframe
Correct_Colnames <- function(df) {
  delete.columns <- grep("(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
  if (length(delete.columns) > 0) {
    row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
    df <- df[,-delete.columns]
    colnames(df) <- gsub("^X", "",  colnames(df))
  }
  return(df)
}

# Factor a list of given columns
factorthese <- function(df, somecolumns){
  Fctr <- names(df) %in% somecolumns
  df[,Fctr] <- lapply(df[,Fctr], function(column) as.factor(as.character(column)))
  return(df)
}

# Transform data (doesn't work that well)
transformit <- function(df){
  df1 <- data.frame(Patient.ID = row.names(df), df, row.names=NULL) 
}

# Create patient ID from sample ID
samptopat <- function(splist){
  splist <- as.character(splist)
  lslist <- vector(mode="character",length = length(splist))
  samples <- regexpr("^TCGA",splist)
  spselect <- samples != -1
  lslist[!spselect] <- splist[!spselect]
  sp_pos_fix <- regexpr("TCGA[[:punct:]]{1}[[:alnum:]]{2}[[:punct:]]{1}[[:alnum:]]{4}",
                        splist[spselect],
                        perl = TRUE)
  lslist[spselect] <- substr(splist[spselect],
                             sp_pos_fix,sp_pos_fix+attributes(sp_pos_fix)[[1]]-1)
  return(factor(lslist))
}

# Splitting Amino acid changes
## Take last letter
lastvalue <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Find position of AA
AApos <- function(splist){
  splist <- as.character(splist)
  lslist <- vector(mode="character",length = length(splist))
  samples <- regexpr("[[:digit:]]",splist)
  spselect <- samples != -1
  lslist[!spselect] <- splist[!spselect]
  sp_pos_fix <- regexpr("[[:digit:]]{1,10}",
                        splist[spselect],
                        perl = TRUE)
  lslist[spselect] <- substr(splist[spselect],
                             sp_pos_fix,sp_pos_fix+attributes(sp_pos_fix)[[1]]-1)
  return(factor(lslist))
}

# Remove NA from columns specified
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Clean Genesets
takegenelevels <- function(df){
  df1 <- tail(df, -1)
  geneset <- as.character(df1[, 1])
  return(geneset)
}

# PCA functions
## Find the correlation between variable and PC
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

## Calculate the contribution of variable to the PC
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}


# Calculate gene means
CalGeneMean <- function(df, genecol){
  genes <-    df[eval(substitute(genecol), df) == T, ]
  genes1 <- droplevels(subset(genes, DataType == "MicroarrayZScore"))
  GeneScore <- dcast(genes1, Hugo_Symbol ~ Patient.ID, value.var = "Value")
  GeneMean <- colMeans(GeneScore[-1], na.rm = T)
  GeneData <- melt(GeneMean)
  GeneData1 <- cbind(rownames(GeneData), data.frame(GeneData, row.names=NULL))
  colnames(GeneData1)[colnames(GeneData1) == "rownames(GeneData)"] <- "Patient.ID"
  GeneData2 <- GeneData1
  return(GeneData2)
}

# Add regression formula
lm_eqn <- function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


# Plots plot
# PC Score Plots
ggComponent <- function(df){
  p <- ggplot(df,aes(x = subgroup, y = Rank)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(subgroup, fill = subgroup),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Cell Type",
         y = "Rank of Component Score"
    )+ 
    geom_dotplot(binaxis='y',
                 method="histodot",
                 stackdir='center',
                 binwidth = 20,
                 position=position_jitter(0.1),
                 alpha=0,
                 dotsize=0.4)+
    theme_bw()+
    theme(
      axis.text = element_text(size=16)) +
    theme(legend.direction = 'horizontal', 
          legend.position = 'top')
  sig <- p+stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif")
  return(sig)}
