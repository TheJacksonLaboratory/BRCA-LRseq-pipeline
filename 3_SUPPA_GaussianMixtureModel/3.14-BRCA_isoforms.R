# List of Differentially spliced subpopulations (panel 4, Fig. 4A)
# Gaussian mixture modeling analysis to find cancer splicing events - TCGA + GTEX
library(rtracklayer)
library(dplyr)

rm(list = ls())


source("3.2-BRCA_isoforms.R")

AL <- read.table(file = "suppa_250bp/AL_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
AF <- read.table(file = "suppa_250bp/AF_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
SE <- read.table(file = "suppa_gtex/SE_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
MX <- read.table(file = "suppa_gtex/MX_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
RI <- read.table(file = "suppa_gtex/RI_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
A5 <- read.table(file = "suppa_gtex/A5_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
A3 <- read.table(file = "suppa_gtex/A3_clusters_annotated.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

all_clusters <- rbind(
  A3,
  A5,
  AF,
  AL,
  SE,
  RI,
  MX)

# Filter duplicate events based on transcript mapping
# Event maps to same transcripts and it is the same type type

all_clusters_ndup <- all_clusters %>% distinct(alternative_transcripts, alternative2_transcripts, type,
                           .keep_all = TRUE)

# Add transcriptome origin to splice events
all_clusters_ndup$transcriptome_origin <- getTranscriptomeOrigin(all_clusters_ndup)
table(all_clusters_ndup$transcriptome_origin)


all_clusters_ndup %>% group_by(type) %>%
  dplyr::summarise(n())


# Add deltaPSI for adjacent TCGA and Gtex

all_clusters_ndup$absDeltaPSI_TvsAdjacent <- 
  abs(all_clusters_ndup$meanPSI_T.x - all_clusters_ndup$mean_PSI_adjacent_TCGA)
all_clusters_ndup$absDeltaPSI_TvsGtex <- 
  abs(all_clusters_ndup$meanPSI_T.x - all_clusters_ndup$mean_PSI_gtex_all)

range(all_clusters_ndup$absDeltaPSI_TvsAdjacent, na.rm = TRUE)
range(all_clusters_ndup$absDeltaPSI_TvsGtex, na.rm = TRUE)
range(all_clusters_ndup$absDeltaPSI_TvsN.x)

#********* P-value adjustment using BH method (a.k.a FDR)
all_clusters_ndup$Wilcox_TvsN_adjusted_p <- 
  p.adjust(all_clusters_ndup$Wilcox_pval.x, method = "BH")

range(all_clusters_ndup$Wilcox_pval.x)
range(all_clusters_ndup$Wilcox_TvsN_adjusted_p)

all_clusters_sig <- all_clusters_ndup %>% filter(Counts.x >= 50,
                                all_clusters_ndup$Wilcox_TvsN_adjusted_p < 0.01,                 
                                absDeltaPSI_TvsAdjacent >= 0.20,
                                absDeltaPSI_TvsGtex >= 0.20)

# all_clusters_sig <- all_clusters_ndup %>% filter(Counts.x >= 50,
#                                                  all_clusters_ndup$Wilcox_TvsN_adjusted_p < 0.01,                 
#                                                  absDeltaPSI_TvsN.x >= 0.2)

all_clusters_sig %>% group_by(type) %>%
  dplyr::summarise(n())

write.table(all_clusters_sig, 
            file = "suppa_gtex/all_clusters_differentialSplicing.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)


all_clusters_sig <- read.table(file = "suppa_gtex/all_clusters_differentialSplicing.tsv", 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(all_clusters_sig)

all_clusters_sig %>% group_by(transcriptome_origin) %>%
  dplyr::summarise(n())
