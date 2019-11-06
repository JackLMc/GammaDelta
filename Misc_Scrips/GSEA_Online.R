load("Bulk/Final_edgeR.RData")

### Write out
output_dir <- "/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/4_Gamma_Delta/GSEA"

brca_hd_tep_tmm_normalized_expression_df <- x$counts

expression_dataset_path <- file.path(output_dir, "x.txt")
write.table(brca_hd_tep_tmm_normalized_expression_df,
            quote=FALSE,
            sep = "\t",
            file=expression_dataset_path,
            row.names = FALSE)

### Calculate variability (dispersions) in data
brca_hd_tep_tmm_normalized_dge <- x
brca_hd_tep_fitted_commondisp_dge <- edgeR::estimateCommonDisp(brca_hd_tep_tmm_normalized_dge)
brca_hd_tep_fitted_tagwise_dge <- edgeR::estimateTagwiseDisp(brca_hd_tep_fitted_commondisp_dge)

### Perform differential expression testing (comparison is 'BrCa' vs 'HD')
baseline_class <- "VD1.CD27HI"
test_class <- "VD1.CD27LO"
comparison <- c(baseline_class, test_class)
brca_hd_tep_de_tested_dge <- edgeR::exactTest(brca_hd_tep_fitted_tagwise_dge, pair = comparison)

### Perform multiple-testing correction using Benjamini-Hockberg procedure
brca_hd_tep_de_tested_tt <- edgeR::topTags(brca_hd_tep_de_tested_dge,
                                           n = nrow(brca_hd_tep_tmm_normalized_dge),
                                           adjust.method = "BH",
                                           sort.by = "PValue")

### Rank by inverse of p-value taking into account 'sign' of change in BrCa (i.e. increase/decrease) relative to HD
brca_hd_tep_rank_values <- sign(brca_hd_tep_de_tested_tt[["table"]][["logFC"]]) * (-1) * log10(brca_hd_tep_de_tested_tt[["table"]][["PValue"]])

### Take into account log10(0) = -Inf
brca_hd_tep_rank_values_max <- max(brca_hd_tep_rank_values[ brca_hd_tep_rank_values != Inf ])
brca_hd_tep_rank_values_unique <- sapply( brca_hd_tep_rank_values,
                                          function(x) replace(x, is.infinite(x),
                                                              sign(x) * (brca_hd_tep_rank_values_max + runif(1))) )


### Construct the data frame we wish place into a tabular file
genenames <- (rownames(brca_hd_tep_de_tested_tt[["table"]]))
brca_hd_tep_ranks_df <- data.frame(gene=genenames,
                                   rank=brca_hd_tep_rank_values_unique,
                                   stringsAsFactors = FALSE)
brca_hd_tep_ordered_ranks_df <- brca_hd_tep_ranks_df[order(brca_hd_tep_ranks_df[,2], decreasing = TRUE), ]

## Write out to file
rank_list_path <- file.path(output_dir, "brca_hd_tep.rnk")
write.table(brca_hd_tep_ordered_ranks_df,
            quote=FALSE,
            sep = "\t",
            file=rank_list_path,
            row.names = FALSE)


brca_hd_tep_ordered_ranks_file_df <- read.table(rank_list_path,
                                                header = TRUE,
                                                check.names = FALSE)

ranks_head <- head(brca_hd_tep_ordered_ranks_file_df, n=5)
rownames(ranks_head) <- NULL
knitr::kable(ranks_head)

brca_hd_tep_filtered_dge <- x
n_samples <- dim(brca_hd_tep_filtered_dge)[2]
n_classes <- 2

l1 <- paste(n_samples, n_classes, "1")
l2 <- paste("#", brca_hd_tep_de_tested_tt[["comparison"]][1], brca_hd_tep_de_tested_tt[["comparison"]][2])
l3 <- paste(brca_hd_tep_filtered_dge[["samples"]][["group"]], collapse = " ")
brca_hd_tep_cls <- rbind(l1, l2, l3)

### Write out to file
categorical_class_path <- file.path(output_dir, "brca_hd_tep.cls")
write(brca_hd_tep_cls,
      file=categorical_class_path,
      sep = "\t")
rownames(brca_hd_tep_cls) <- NULL
brca_hd_tep_cls



doEnrichment <- T
### Declare user-defined settings
gsea_jar_path <- file.path("/Users/jlm650/Downloads/gsea-3.0.jar")
gsea_rpt_label <- "tep_BrCa_HD_analysis"
gsea_analysis_name <- "tep_BrCa_HD"
gsea_out <- file.path(getwd(), "gsea_output")
gsea_gmx <- file.path(getwd(),
                      "data",
                      "Human_GOBP_AllPathways_no_GO_iea_August_01_2017_symbol.gmt")
gsea_rank_list_path <- rank_list_path
gsea_num_permutations <- 1000
gsea_min_gs_size <- 15
gsea_max_gs_size <- 200

## Execute GSEA
command <- paste("java -cp", gsea_jar_path,
                 "-Xmx1G xtools.gsea.GseaPreranked",
                 "-rpt_label", gsea_analysis_name,
                 "-out", gsea_out,
                 "-gmx", gsea_gmx,
                 "-rnk", gsea_rank_list_path,
                 "-nperm", gsea_num_permutations,
                 "-set_min", gsea_min_gs_size,
                 "-set_max", gsea_max_gs_size,
                 "-collapse false",
                 "-scoring_scheme weighted",
                 "-permute gene_set",
                 "-num 100",
                 "-plot_top_x 20",
                 "-rnd_seed 12345",
                 "-zip_report false",
                 "-gui false",
                 ">", paste("gsea_output_", gsea_rpt_label, ".txt", sep=""),
                 sep=" ")

if( isTRUE ( doEnrichment) ){
  system(command)
}


gsea_tep_BrCa_HD_analysis_directory <- list.files(path = gsea_out, pattern = "\\.GseaPreranked")
