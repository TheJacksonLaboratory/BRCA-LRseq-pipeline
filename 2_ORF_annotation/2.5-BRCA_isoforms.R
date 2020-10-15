library(dplyr)
library(ggplot2)
library(ggsci)
library(scales)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(parallel)

# Transdecoder precision analysis on FSM transcripts with exact first and last exons
# ***** RUN on cluster using 2.5-BRCA_isoforms.pbs

# TP - fsm exact ORFs > 99% identity
# FP - fsm exact ORFs < 99% identity
# Precision = (TP) / (TP+FP)

pac_gtf <- import.gff("/projects/dveiga/analysis/git/BRCA_isoforms/QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf.gz")
genc_gtf <- import.gff("/projects/banchereau-lab/tools/Gencode_30/gencode.v30.annotation.gtf")

pac_gtf <- import.gff("QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf.gz")
genc_gtf <- import.gff("/Users/dveiga/tools/Gencode/gencode.v30.annotation.sorted.gff.gz")


exons_genc <- genc_gtf %>% as.data.frame() %>%
  filter(type %in% "exon") 
exons_pac <- pac_gtf %>% as.data.frame() %>%
  filter(type %in% "exon")

# Sqanti transcript annotation QC pass
sq_QCpass <- read.table("./QC_passed/Sqanti_annotation_QC_pass_transcripts_unique.txt",
                        sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)

sq_QCpass_fsm <- sq_QCpass %>% filter(structural_category %in% "full-splice_match")
head(sq_QCpass_fsm)

res <- mclapply(seq(1:nrow(sq_QCpass_fsm)), function(i){
  
  exons_genc_txn <- exons_genc %>% 
    filter(transcript_id %in% sq_QCpass_fsm$associated_transcript[i]) 
  exons_genc_txn <- as(exons_genc_txn, "GRanges")
  
  exons_pac_txn <- exons_pac %>%
    filter(transcript_id %in% sq_QCpass_fsm$isoform[i])
  exons_pac_txn <- as(exons_pac_txn, "GRanges")
  
  ov_first <- findOverlaps(exons_genc_txn[1], exons_pac_txn[1], 
                           type = "equal", minoverlap = 1, maxgap = 10)
  
  ov_last <- findOverlaps(exons_genc_txn[length(exons_genc_txn)], exons_pac_txn[length(exons_pac_txn)],
                          type = "equal", minoverlap = 1, maxgap = 10)
  

  return(data.frame(isoform = sq_QCpass_fsm$isoform[i],
                    exact_first =   length(queryHits(ov_first)),
                    exact_last =   length(queryHits(ov_last))
                      ))

}, mc.cores = 6)
  
  
df_res <- do.call(rbind, res)
save(df_res, file = "./transdecoder_precision/results_first_last_Exact.Rd")
head(df_res)

df_select <- df_res %>% filter(exact_first %in% 1, exact_last %in% 1)

df_res <- df_res %>% filter(isoform %in% sq_QCpass$isoform)

# ************ Compute precision

# Approach 1 - use only FSM transcripts for which First/exon match the reference
#load("./transdecoder_precision/results_first_last_Exact.Rd")

df_select <- df_res %>% filter(exact_first %in% 1, exact_last %in% 1)
head(df_select)

load("./transdecoder_globalAlignment/Transdecoder_map_aligned.Rd")
#load("./transdecoder_globalAlignment/Transdecoder_map_aligned_PID4.Rd")
transdecoder_aligned <- transdecoder_aligned %>% filter(PBid %in% sq_QCpass$isoform)
head(transdecoder_aligned)


df_select_annot <- inner_join(df_select, transdecoder_aligned, by = c("isoform" = "PBid" ))

head(df_select_annot)

TP <- df_select_annot %>% filter(p_identity == 100,
                                 Transdecoder_ORF_type %in% "complete") %>% dplyr::count()
FP <- df_select_annot %>% filter(p_identity < 100,
                                 Transdecoder_ORF_type %in% "complete") %>% dplyr::count()

prec = TP / (TP+FP)
prec

# Approach 2 - Using all FSM transcripts
df_select_annot <- inner_join(df_res, transdecoder_aligned, by = c("isoform" = "PBid" ))
TP <- df_select_annot %>% filter(p_identity > 99.9,
                                 Transdecoder_ORF_type %in% "complete") %>% dplyr::count()
FP <- df_select_annot %>% filter(p_identity <= 99.9,
                                 Transdecoder_ORF_type %in% "complete") %>% dplyr::count()

prec = TP / (TP+FP)
prec


