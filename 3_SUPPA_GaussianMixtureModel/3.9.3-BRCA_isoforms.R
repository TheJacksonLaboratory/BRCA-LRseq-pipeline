# Compute Survival within clusters - SE, RI, MX, A3, A5, AL, AF events
# Using GMM subpopulations

rm(list = ls())

source("3.2-BRCA_isoforms.R")

library(rtracklayer)
library(dplyr)


# ********* A5 events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A5_clust_BRCA_TCGA.RData")
dim(A5_clust$ts_clusters)

out_folder <- "BRCA_TCGA_A5_subpopulation"

head(A5_clust$ts_clusters)
# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis_by_subpop(A5_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 A5_clust$ts_clusters, out_folder)

# Significant events: Global p-value < 0.01 and adjusted pairwise comparison < 0.05
surv_df %>% filter(Survival_overall_sig == TRUE) %>% nrow()

# # Add survival p-value to clust object
A5_ts_clusters_surv <- inner_join(A5_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
A5_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>% count()
A5_ts_clusters_surv$type <- "A5"


# ********* A3 events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A3_clust_BRCA_TCGA.RData")
dim(A3_clust$ts_clusters)

out_folder <- "BRCA_TCGA_A3_subpopulation"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis_by_subpop(A3_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 A3_clust$ts_clusters, out_folder)

# Significant events: Global p-value < 0.01 and adjusted pairwise comparison < 0.05
surv_df %>% filter(Survival_overall_sig == TRUE) %>% nrow()

# # Add survival p-value to clust object
A3_ts_clusters_surv <- inner_join(A3_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
A3_ts_clusters_surv$type <- "A3"

# ********* MX events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/MX_clust_BRCA_TCGA.RData")
dim(MX_clust$ts_clusters)

out_folder <- "BRCA_TCGA_MX_subpopulation"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis_by_subpop(MX_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 MX_clust$ts_clusters, out_folder)

# # Add survival p-value to clust object
MX_ts_clusters_surv <- inner_join(MX_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
MX_ts_clusters_surv$type <- "MX"

# ********* SE events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/SE_clust_BRCA_TCGA.RData")
dim(clust$ts_clusters)

out_folder <- "BRCA_TCGA_SE_subpopulation"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis_by_subpop(clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 clust$ts_clusters, out_folder)

# # Add survival p-value to clust object
SE_ts_clusters_surv <- inner_join(clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
SE_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>% count()
SE_ts_clusters_surv$type <- "SE"

# ********* RI events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/RI_clust_BRCA_TCGA.RData")
dim(RI_clust$ts_clusters)

out_folder <- "BRCA_TCGA_RI_subpopulation"

# Compute survival p-values for ts_clusters
surv_df <- BRCA_survivalAnalysis_by_subpop(RI_clust,
                                 clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                 RI_clust$ts_clusters, out_folder)

# # Add survival p-value to clust object
RI_ts_clusters_surv <- inner_join(RI_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
RI_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>% count()
RI_ts_clusters_surv$type <- "RI"


# # ********* AF events
# # Load clust object with GMM results
# load("./suppa_gtex/gmm_results/AF_clust_BRCA_TCGA.RData")
# dim(AF_clust$ts_clusters)
# 
# out_folder <- "BRCA_TCGA_AF"
# 
# # Compute survival p-values for ts_clusters
# surv_df <- BRCA_survivalAnalysis_by_subpop(AF_clust,
#                                  clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
#                                  AF_clust$ts_clusters, out_folder)
# 
# # # Add survival p-value to clust object
# AF_ts_clusters_surv <- inner_join(AF_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
# AF_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>% count()
# AF_ts_clusters_surv %>% filter(is.na(Survival_pval)) %>% count()
# AF_ts_clusters_surv$type <- "AF"
# 
# # ********* AL events
# # Load clust object with GMM results
# load("./suppa_gtex/gmm_results/AL_clust_BRCA_TCGA.RData")
# dim(AL_clust$ts_clusters)
# 
# out_folder <- "BRCA_TCGA_AL"
# 
# # Compute survival p-values for ts_clusters
# surv_df <- BRCA_survivalAnalysis_by_subpop(AL_clust,
#                                  clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
#                                  AL_clust$ts_clusters, out_folder)
# 
# # # Add survival p-value to clust object
# AL_ts_clusters_surv <- inner_join(AL_clust$ts_clusters, surv_df, by = c("Junction", "Cluster", "Gene") )
# AL_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>% count()
# AL_ts_clusters_surv %>% filter(is.na(Survival_pval)) %>% count()
# AL_ts_clusters_surv$type <- "AL"
# 
# head(A5_ts_clusters_surv)

all_ts_clusters_surv <- rbind(A5_ts_clusters_surv, A3_ts_clusters_surv,
                              MX_ts_clusters_surv, SE_ts_clusters_surv,
                              RI_ts_clusters_surv)
dim(all_ts_clusters_surv)

all_ts_clusters_surv %>% filter(Survival_overall_sig == TRUE) %>%
  group_by(type) %>% count()

all_ts_clusters_surv %>% group_by(type) %>% count()

save(all_ts_clusters_surv, file = "suppa_gtex/all_events_ts_clusters_surv_subpopulation.RData")




