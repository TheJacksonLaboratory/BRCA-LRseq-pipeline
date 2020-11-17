# Compute Survival within clusters -AL, AF events 250bp
# Using GMM subpopulations

rm(list = ls())

source("3.2-BRCA_isoforms.R")

library(rtracklayer)
library(dplyr)

# ********* AF events
# Load clust object with GMM results
load("./suppa_250bp/AF_250bp_clust_BRCA_TCGA.RData")
dim(AF_250bp_clust$ts_clusters)

out_folder <- "BRCA_TCGA_AF_250bp"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis(AF_250bp_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 AF_250bp_clust$ts_clusters, out_folder)

# # Add survival p-value to clust object
AF_250bp_ts_clusters_surv <- inner_join(AF_250bp_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
AF_250bp_ts_clusters_surv %>% filter(Survival_pval < 0.01) %>% count()
AF_250bp_ts_clusters_surv %>% filter(is.na(Survival_pval)) %>% count()
AF_250bp_ts_clusters_surv$type <- "AF_250bp"

# ********* AL events
# Load clust object with GMM results
load("./suppa_250bp/AL_250bp_clust_BRCA_TCGA.RData")
dim(AL_250bp_clust$ts_clusters)

out_folder <- "BRCA_TCGA_AL_250bp"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis(AL_250bp_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 AL_250bp_clust$ts_clusters, out_folder)

# # Add survival p-value to clust object
AL_250bp_ts_clusters_surv <- inner_join(AL_250bp_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
AL_250bp_ts_clusters_surv %>% filter(Survival_pval < 0.01) %>% count()
AL_250bp_ts_clusters_surv %>% filter(is.na(Survival_pval)) %>% count()
AL_250bp_ts_clusters_surv$type <- "AL_250bp"


all_ts_clusters_surv_250 <- rbind(AF_250bp_ts_clusters_surv,
                              AL_250bp_ts_clusters_surv)
dim(all_ts_clusters_surv_250)

all_ts_clusters_surv_250 %>% filter(Survival_pval < 0.01) %>%
  group_by(type) %>% count()

all_ts_clusters_surv_250 %>% group_by(type) %>% count()

save(all_ts_clusters_surv_250, file = "suppa_250bp/all_events_ts_clusters_surv_250bp.RData")
