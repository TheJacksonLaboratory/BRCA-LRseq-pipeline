
# Gaussian mixture modeling analysis to find cancer splicing events
# TCGA + GTEX
# ********** AF events ****************

source("3.2-BRCA_isoforms.R")

library(rtracklayer)
library(dplyr)

# # ********* Cluster variables
gff_file <- "/projects/dveiga/analysis/git/BRCA_isoforms/suppa/PacBio_NIC_NNC_ISM_Gencode30_merged.gtf"
PSI_file <- "/projects/dveiga/analysis/git/BRCA_isoforms/suppa/AF_PSI.psi"
PSI_file_gtex <- "/projects/dveiga/analysis/git/BRCA_isoforms/suppa_gtex/AF_PSI_gtex.psi"
events_file <- "/projects/dveiga/analysis/git/BRCA_isoforms/suppa/suppa_events/AS_events_AF_variable_10.ioe"
analysis_name <- "BRCA_TCGA_GTEX"
clust_dir <- "AF_ts3_clust_output/" # Output folder (Rd file with df.exp for each cluster will be saved there)
output_file <- "/projects/dveiga/analysis/git/BRCA_isoforms/suppa_gtex/AF_clust_BRCA_TCGA.RData"


gtf <- import.gff(gff_file)
SEpsi <- import_suppa_TCGA(fpath = PSI_file)
head(SEpsi[1, 1:10])

events <- read.table(file = events_file, sep = "\t",
                     header = TRUE, stringsAsFactors = FALSE)
idx.gene <- match(events$gene_id, gtf$gene_id)
events$gene_symbol <- gtf$gene_name[idx.gene]
head(events)

barcode_annot <- annotate_TCGA_barcode(colnames(SEpsi))

head(barcode_annot)
table(barcode_annot$primarySite)
table(barcode_annot$sample_shortLetterCode)


barcode_tumor <- barcode_annot %>%
  filter(Study_abbrev %in% c("BRCA"),
         sample_shortLetterCode %in% c("TP"))

barcode_control <- barcode_annot %>%
  filter(sample_shortLetterCode %in% c("NT"))

PSI_tumor <- SEpsi[ , colnames(SEpsi) %in% barcode_tumor$Barcode]
PSI_control <- SEpsi[ , colnames(SEpsi) %in% barcode_control$Barcode]

# Sort PSI_tumor and PSI_control using order from events
idx.events <- match(events$event_id, rownames(PSI_tumor))
sum(is.na(idx.events))
PSI_tumor <- PSI_tumor[idx.events, ]
PSI_control <- PSI_control[idx.events, ]

head(rownames(PSI_tumor))
head(events$event_id)


# Load GTEx controls
SEpsi_gtex <- import_suppa_GTEX(fpath = PSI_file_gtex)
head(SEpsi_gtex[1, 1:10])

gtex_annot <- annotate_GTEX_barcode(colnames(SEpsi_gtex))

table(gtex_annot$histological_type)

# Sort PSI_tumor and PSI_control using order from events
idx.events.gtex <- match(events$event_id, rownames(SEpsi_gtex))
sum(is.na(idx.events.gtex))
PSI_gtex <- SEpsi_gtex[idx.events.gtex, ]

head(rownames(PSI_gtex))
head(rownames(PSI_control))
head(events$event_id)


#Parameters for clustering
min_cluster = 10 #clusters below this threshold are merged
min_tumor = 30 #min number of tumor samples to be clustered after removing NA samples
min_controls = 10 #min number of control samples to be clustered after removing NA

# 
# test <- 1000
# PSI_tumor <- PSI_tumor[1:test, ]
# PSI_control <- PSI_control[1:test, ]
# PSI_gtex <- PSI_gtex[1:test, ]
# events <- events[1:test, ]

PSI_all_control <- cbind(PSI_control, PSI_gtex)

mixtureClusterAnalysis(analysis_name, PSI_tumor, PSI_all_control,
                       events, min_cluster, min_tumor, min_controls,
                       clust_dir)

# nonParallel_mixtureClusterAnalysis(analysis_name, PSI_tumor, PSI_control,
#                                    events, min_cluster, min_tumor, min_controls,
#                                    clust_dir)

#Collect results 
clust <- list()
clust$composition <- computeClusterComposition(analysis_name, clust_dir, events)
clust$meanPSI <- computeClusterMeanPSI(analysis_name, clust_dir, events)

params <- list(
  min_size = 30, #number of patients
  tumor_purity = 90, # percent of tumor purity
  deltaPSI = 0.1, #minimum deltaPSI
  pvalue = 0.01 #minimum pvalue  
)

clust$ts_clusters <- selectClusters(clust, params, clust_dir, analysis_name)
length(unique(clust$ts_clusters$Gene)) #number of unique genes 
length(unique(clust$ts_clusters$Junction)) #number of exon skipping events

#Recover patient assigment and expression for TS clusters 
clust$cluster_exp <- gatherClusterExp(analysis_name, clust_dir, events, 
                                      clust$ts_clusters)

AF_clust <- clust
#Save results
save(AF_clust, file = output_file)

