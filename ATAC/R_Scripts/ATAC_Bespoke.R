# library(Rsubread)
# 
# genome <- "./Bulk/AME/hg19/uncompressed/hg19_ref_genome.fa"
# indexForSubread <- gsub("\\.fa$", "", genome)
# 
# buildindex(indexForSubread, genome, indexSplit = FALSE)

# Look at proportion of mapped reads.
for(i in length(files)){
  sortedBAM <- i
  pmapped <- propmapped(sortedBAM)
  this_list[[i]] <- pmapped
}
sortedBAM <- "/Volumes/ResearchData/Willcox Group/Jack/GD_Transcriptomic_paper/GD_ATAC/bams/CD8-Naive.bam"
library(Rsubread)


pmapped <- propmapped(sortedBAM)
pmapped


library(Rsamtools)
library(ggplot2)
library(magrittr)

idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
  geom_bar(stat = "identity") + coord_flip()



library(GenomicAlignments)
atacReads <- readGAlignmentPairs(sortedBAM, 
                                 param = ScanBamParam(mapqFilter = 1,
                                                      flag = scanBamFlag(isPaired = T, isProperPair = T), 
                                                                 what = c("qname", "mapq", "isize")))



atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)



library(magrittr)
library(dplyr)
library(ggplot2)
fragLenPlot <- table(insertSizes) %>% 
  data.frame %>% 
  rename(InsertSize = insertSizes, Count = Freq) %>%
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>%
  ggplot(aes(x = InsertSize, y = Count)) + 
  geom_line()

fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()

fragLenPlot + scale_y_continuous(trans = "log2") +
  geom_vline(xintercept = c(180, 247), colour = "red") +
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), fix = "start", 1)

# BiocManager::install("soGGi")
library(soGGi)

TSSs <- trim(TSSs)

# Nucleosome free
nucFree <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                      format = "bam", paired = T, minFragmentLength = 0, maxFragmentLength = 100, 
                      forceFragment = 50)

# Mononucleosome
monoNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
                      forceFragment = 80)

# Dinucleosome
diNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                    format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
                    forceFragment = 160)


# BiocManager::install(c("ChIPQC", "rtracklayer", "DT"))
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)


getwd()
blkList <- import.bed("./ATAC/Data/Blcklist.bed.gz")
openRegionPeaks <- "/Volumes/JackMcMurray/ATACSeq/3_peak_called/VD1-Eff_peaks.narrowPeak"

qcRes <- ChIPQCsample("/Volumes/JackMcMurray/ATACSeq/2_bam/VD1-Eff.bam", 
                      peaks = openRegionPeaks, annotation = "hg19", chromosomes = NULL, 
                      blacklist = blkList, verboseT = T)

?ChIPQCsample

QCmetrics(qcRes) %>% t %>%
  data.frame %>% 
  dplyr:::select(Reads, starts_with(c("Filt")),
                 starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)

flagtagcounts(qcRes) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
                                                       100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)





getwd()
peaks <- read.delim("./ATAC/Data/Coverage_allPeaks_400.txt")

levels(peaks$Chr)

meaningful <- droplevels(subset(peaks, Chr == "chr1" | Chr == "chr2" | 
                    Chr == "chr3" | Chr == "chr4" |
                    Chr == "chr5" | Chr == "chr6" |
                    Chr == "chr7" | Chr == "chr8" |
                    Chr == "chr9" | Chr == "chr10" |
                    Chr == "chr11" | Chr == "chr12" |
                    Chr == "chr13" | Chr == "chr14" |
                    Chr == "chr15" | Chr == "chr16" |
                    Chr == "chr17" | Chr == "chr18" |
                    Chr == "chr19" | Chr == "chr20" |
                    Chr == "chr21" | Chr == "chr22" |
                    Chr == "chrX" | Chr == "chrY" | Chr == "chrM"))

library(tidyverse) 
peaks1 <- meaningful %>% gather(contains("_treat_pileup.bdg.bedGraph.avg.over.400.bp"), key = "Population", value = "Peak")
peaks1$Population <- as.factor(peaks1$Population)


this <- droplevels(subset(peaks1, Population == "CD8.CM_treat_pileup.bdg.bedGraph.avg.over.400.bp"))

idxstatsBam(this) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
  geom_bar(stat = "identity") + coord_flip()


head(this)

??read_idxstats

ChIPseeker:: plotAnnoPie(this)

??plotAnnoPie

library(ChIPseeker)
MacsCalls_chr20_filteredAnno <- annotatePeak(MacsCalls_chr20_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)

