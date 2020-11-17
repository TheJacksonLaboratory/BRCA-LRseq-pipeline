
# Master table with new survival analysis
# Annotation of all events with transcripts, tumor ids, control ids, survival, pam50
# Gaussian mixture modeling analysis to find cancer splicing events - TCGA + GTEX
library(rtracklayer)
library(dplyr)
library(stringr)
library(tibble)

rm(list = ls())

source("3.2-BRCA_isoforms.R")

# Load PAM50 enrichment results
load("suppa_gtex/all_events_ts_clusters_pam50.RData")
#save(all_ts_clusters_pam50, file = "suppa_gtex/all_events_ts_clusters_pam50.RData")
pam50 <- all_ts_clusters_pam50 %>% 
  dplyr::select(c("Junction", "Gene", "Cluster",
                  "Basal_Like", "LuminalA", "LuminalB", "HER2")) %>%
  dplyr::rename("Basal_Like_log10_pval" = "Basal_Like",
                "LuminalA_log10_pval" = "LuminalA", 
                "LuminalB_log10_pval" = "LuminalB", 
                "HER2_log10_pval" = "HER2")


#********* P-value adjustment using BH method (a.k.a FDR)
#Basal-like 
adjusted_p <- p.adjust(10^(pam50$Basal_Like_log10_pval*-1), method = "BH" ) 
pam50$Basal_Like_log10_pval_adj <- -1*log10(adjusted_p)
head(pam50$Basal_Like_log10_pval)
head(pam50$Basal_Like_log10_pval_adj)

#LuminalA
adjusted_p <- p.adjust(10^(pam50$LuminalA_log10_pval*-1), method = "BH" ) 
pam50$LuminalA_log10_pval_adj <- -1*log10(adjusted_p)

#LuminalB
adjusted_p <- p.adjust(10^(pam50$LuminalB_log10_pval*-1), method = "BH" ) 
pam50$LuminalB_log10_pval_adj <- -1*log10(adjusted_p)

#HER2
adjusted_p <- p.adjust(10^(pam50$HER2_log10_pval*-1), method = "BH" ) 
pam50$HER2_log10_pval_adj <- -1*log10(adjusted_p)

# Select only Adjusted p-values columns
pam50_adj <- pam50 %>% dplyr::select(c("Junction", "Gene", "Cluster",
                                       "Basal_Like_log10_pval_adj", "LuminalA_log10_pval_adj",
                                       "LuminalB_log10_pval_adj", "HER2_log10_pval_adj"))

head(pam50_adj)

head(all_ts_clusters_pam50)

# Load survival enrichment results
#save(all_ts_clusters_surv, file = "suppa_gtex/all_events_ts_clusters_surv.RData")
load("suppa_gtex/all_events_ts_clusters_surv_subpopulation.RData")
table(all_ts_clusters_surv$type)
head(all_ts_clusters_surv)
# surv <- all_ts_clusters_surv %>% 
#   dplyr::select(c("Junction", "Gene", "Cluster", "type",
#                   "Survival_pval"))
surv <- all_ts_clusters_surv
head(surv)

load("suppa_250bp/all_events_ts_clusters_surv_250bp_subpopulation.RData")
table(all_ts_clusters_surv_250$type)
# surv_250 <- all_ts_clusters_surv_250 %>% 
#   dplyr::select(c("Junction", "Gene", "Cluster", "type",
#                   "Survival_pval"))
surv_250 <- all_ts_clusters_surv_250

# Still contain similar AS events to be removed 
surv %>% filter(Survival_overall_sig == TRUE) %>% 
  dplyr::group_by(type) %>%
  dplyr::summarise(n())

# surv$Survival_pval_adj <- p.adjust(surv$Survival_pval, method = "BH")
# surv %>% filter(Survival_pval_adj < 0.05) %>% 
#   dplyr::group_by(type) %>%
#   dplyr::summarise(n())
# surv_250$Survival_pval_adj <- p.adjust(surv_250$Survival_pval, method = "BH")
# surv <- surv %>% dplyr::select(c("Junction", "Gene", "Cluster", "Survival_pval_adj"))
# surv_250 <- surv_250 %>% dplyr::select(c("Junction", "Gene", "Cluster", "Survival_pval_adj"))


# ********* A5 events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A5_clust_BRCA_TCGA.RData")
dim(A5_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa/suppa_events/AS_events_A5_variable_10.ioe") 
head(ioe)

A5_clusters_annot <- inner_join(A5_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(A5_clust, A5_clusters_annot)

A5_clusters_annot2 <- inner_join(A5_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and GTEx breast)
# Add mean PSI and number of control samples detected per tissue
control_ids <- getControlsInClusters(A5_clust, A5_clusters_annot)

A5_clusters_annot3 <- inner_join(A5_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))


# 4. Add PAM50 annotation
A5_clusters_annot4 <- inner_join(A5_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
A5_clusters_annot5 <- inner_join(A5_clusters_annot4, surv,
                                 by = c("Junction", "Cluster", "Gene"))

write.table(A5_clusters_annot5, 
            file = "suppa_gtex/A5_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)


# ********* A3events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/A3_clust_BRCA_TCGA.RData")
dim(A3_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa/suppa_events/AS_events_A3_variable_10.ioe") 
head(ioe)

A3_clusters_annot <- inner_join(A3_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(A3_clust, A3_clusters_annot)

A3_clusters_annot2 <- inner_join(A3_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(A3_clust, A3_clusters_annot)

A3_clusters_annot3 <- inner_join(A3_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))


# 4. Add PAM50 annotation
A3_clusters_annot4 <- inner_join(A3_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
A3_clusters_annot5 <- inner_join(A3_clusters_annot4, surv,
                                 by = c("Junction", "Cluster", "Gene"))

write.table(A3_clusters_annot5, 
            file = "suppa_gtex/A3_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)


# ********* MX events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/MX_clust_BRCA_TCGA.RData")
dim(MX_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa/suppa_events/AS_events_MX_variable_10.ioe") 
head(ioe)

MX_clusters_annot <- inner_join(MX_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(MX_clust, MX_clusters_annot)

MX_clusters_annot2 <- inner_join(MX_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(MX_clust, MX_clusters_annot)

MX_clusters_annot3 <- inner_join(MX_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))


# 4. Add PAM50 annotation
MX_clusters_annot4 <- inner_join(MX_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
MX_clusters_annot5 <- inner_join(MX_clusters_annot4, surv,
                                 by = c("Junction", "Cluster", "Gene"))

write.table(MX_clusters_annot5, 
            file = "suppa_gtex/MX_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)


# ********* RI events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/RI_clust_BRCA_TCGA.RData")
dim(RI_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa/suppa_events/AS_events_RI_variable_10.ioe") 
head(ioe)

RI_clusters_annot <- inner_join(RI_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(RI_clust, RI_clusters_annot)

RI_clusters_annot2 <- inner_join(RI_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(RI_clust, RI_clusters_annot)

RI_clusters_annot3 <- inner_join(RI_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))


# 4. Add PAM50 annotation
RI_clusters_annot4 <- inner_join(RI_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
RI_clusters_annot5 <- inner_join(RI_clusters_annot4, surv,
                                 by = c("Junction", "Cluster", "Gene"))

write.table(RI_clusters_annot5, 
            file = "suppa_gtex/RI_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)

# ********* SE events
# Load clust object with GMM results
load("./suppa_gtex/gmm_results/SE_clust_BRCA_TCGA.RData")
dim(clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa/suppa_events/AS_events_SE_variable_10.ioe") 
head(ioe)

SE_clusters_annot <- inner_join(clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(clust, SE_clusters_annot)

SE_clusters_annot2 <- inner_join(SE_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(clust, SE_clusters_annot)

SE_clusters_annot3 <- inner_join(SE_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))


# 4. Add PAM50 annotation
SE_clusters_annot4 <- inner_join(SE_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
SE_clusters_annot5 <- inner_join(SE_clusters_annot4, surv,
                                 by = c("Junction", "Cluster", "Gene"))

write.table(SE_clusters_annot5, 
            file = "suppa_gtex/SE_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)

# ********* AL events
# Load clust object with GMM results
load("./suppa_250bp/AL_250bp_clust_BRCA_TCGA.RData")
dim(AL_250bp_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa_250bp/suppa_events/AS_events_AL_variable_250.ioe") 
head(ioe)

AL_clusters_annot <- inner_join(AL_250bp_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(AL_250bp_clust, AL_clusters_annot)

AL_clusters_annot2 <- inner_join(AL_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(AL_250bp_clust, AL_clusters_annot)

AL_clusters_annot3 <- inner_join(AL_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))

head(control_ids)

# 4. Add PAM50 annotation
AL_clusters_annot4 <- inner_join(AL_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
AL_clusters_annot5 <- inner_join(AL_clusters_annot4, surv_250,
                                 by = c("Junction", "Cluster", "Gene"))


head(AL_clusters_annot5)

write.table(AL_clusters_annot5, 
            file = "suppa_250bp/AL_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)


# ********* AL events
# Load clust object with GMM results
load("./suppa_250bp/AF_250bp_clust_BRCA_TCGA.RData")
dim(AF_250bp_clust$ts_clusters)

# 1. Add transcript annotation to clusters
ioe <- getSuppaTranscripts("./suppa_250bp/suppa_events/AS_events_AF_variable_250.ioe") 
head(ioe)

AF_clusters_annot <- inner_join(AF_250bp_clust$ts_clusters, ioe, 
                                by = c("Junction" = "event_id") ) %>% distinct()

# 2. Add tumor ids to clusters
tumor_ids <- getPatientInClusters(AF_250bp_clust, AF_clusters_annot)

AF_clusters_annot2 <- inner_join(AF_clusters_annot, tumor_ids,
                                 by = c("Junction", "Cluster", "Gene"))

# 3. Add control IDs and control PSI to cluster annotation (Adjacent TCGA and Gtex breast)
control_ids <- getControlsInClusters(AF_250bp_clust, AF_clusters_annot)

AF_clusters_annot3 <- inner_join(AF_clusters_annot2, control_ids,
                                 by = c("Junction", "Cluster", "Gene"))

head(control_ids)

# 4. Add PAM50 annotation
AF_clusters_annot4 <- inner_join(AF_clusters_annot3, pam50_adj,
                                 by = c("Junction", "Cluster", "Gene"))

# 5. Add survival annotation
AF_clusters_annot5 <- inner_join(AF_clusters_annot4, surv_250,
                                 by = c("Junction", "Cluster", "Gene"))


head(AL_clusters_annot5)

write.table(AF_clusters_annot5, 
            file = "suppa_250bp/AF_clusters_annotated.tsv", sep = "\t", quote = TRUE,
            row.names = FALSE)

