# PAM50 Enrichment within clusters - SE, RI, MX, A3, A5, AL, AF events
# Gaussian mixture modeling analysis to find cancer splicing events - TCGA + GTEX

rm(list = ls())

source("3.2-BRCA_isoforms.R")

library(dplyr)

# ********* A5 events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A5_clust_BRCA_TCGA.RData")
dim(A5_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(A5_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          A5_clust$ts_clusters)

A5_ts_clusters_pam50 <- inner_join(A5_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
A5_ts_clusters_pam50$type <- "A5"
A5_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* A3 events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A3_clust_BRCA_TCGA.RData")
dim(A3_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(A3_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          A3_clust$ts_clusters)

A3_ts_clusters_pam50 <- inner_join(A3_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
A3_ts_clusters_pam50$type <- "A3"
A3_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* MX events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/MX_clust_BRCA_TCGA.RData")
dim(MX_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(MX_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          MX_clust$ts_clusters)

MX_ts_clusters_pam50 <- inner_join(MX_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
MX_ts_clusters_pam50$type <- "MX"
MX_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* SE events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/SE_clust_BRCA_TCGA.RData")
dim(clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          clust$ts_clusters)

SE_ts_clusters_pam50 <- inner_join(clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
SE_ts_clusters_pam50$type <- "SE"
SE_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* RI events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/RI_clust_BRCA_TCGA.RData")
dim(RI_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(RI_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          RI_clust$ts_clusters)

RI_ts_clusters_pam50 <- inner_join(RI_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
RI_ts_clusters_pam50$type <- "RI"
RI_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* AF events 250bp
# Load clust object with GMM results
load("./suppa_250bp/AF_250bp_clust_BRCA_TCGA.RData")
dim(AF_250bp_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(AF_250bp_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          AF_250bp_clust$ts_clusters)

AF_250bp_ts_clusters_pam50 <- inner_join(AF_250bp_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
AF_250bp_ts_clusters_pam50$type <- "AF_250bp"
AF_250bp_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

# ********* AF events 250bp
# Load clust object with GMM results
load("./suppa_250bp/AL_250bp_clust_BRCA_TCGA.RData")
dim(AL_250bp_clust$ts_clusters)

# PAM50 enrichment analysis
pam50 <- enrichment_pam50(AL_250bp_clust, 
                          clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                          AL_250bp_clust$ts_clusters)

AL_250bp_ts_clusters_pam50 <- inner_join(AL_250bp_clust$ts_clusters, pam50, by = c("Junction", "Cluster", "Gene") )
AL_250bp_ts_clusters_pam50$type <- "AL_250bp"
AL_250bp_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>% count()

all_ts_clusters_pam50 <- rbind(A5_ts_clusters_pam50, A3_ts_clusters_pam50,
                               MX_ts_clusters_pam50, SE_ts_clusters_pam50,
                               RI_ts_clusters_pam50, AF_250bp_ts_clusters_pam50,
                               AL_250bp_ts_clusters_pam50)

dim(all_ts_clusters_pam50)

# Basal-like enrichment
all_ts_clusters_pam50 %>% filter(Basal_Like > -1*log10(0.01)) %>%
  group_by(type) %>% count()

# Luminal A enrichment
all_ts_clusters_pam50 %>% filter(LuminalA > -1*log10(0.01)) %>%
  group_by(type) %>% count()

# Luminal B enrichment
all_ts_clusters_pam50 %>% filter(LuminalB > -1*log10(0.01)) %>%
  group_by(type) %>% count()

# Her2 enrichment
all_ts_clusters_pam50 %>% filter(HER2 > -1*log10(0.01)) %>%
  group_by(type) %>% count()

all_ts_clusters_pam50 %>% group_by(type) %>% count()

save(all_ts_clusters_pam50, file = "suppa_gtex/all_events_ts_clusters_pam50.RData")

