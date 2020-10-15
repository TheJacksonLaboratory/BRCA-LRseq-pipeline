library(dplyr)

# Statistics on the Subset of non-NMD and novel ORFs used for functional analysis (PFAM, transmembrane, localization)


load("./transdecoder_globalAlignment/Transdecoder_map_aligned.Rd")
head(transdecoder_aligned)

# NMD sensitive prediction 
load("./ORF_effect/df_PTC_effect.Rd")
head(df_PTC_effect)
non_NMD <- df_PTC_effect %>% filter(PTC_status == "no")

# ********** Filter QC passing non-redundant transcripts
sq_QCpass <- read.table("./QC_passed/Sqanti_annotation_QC_pass_transcripts_unique.txt",
                        sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)

transdecoder_aligned <- transdecoder_aligned %>% filter(PBid %in% sq_QCpass$isoform)

total <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity >= 0
  ) %>% dplyr::count()

novel <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity >= 0,
         p_identity <= 99,
  ) %>% dplyr::count()

known <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity > 99
  ) %>% dplyr::count()

novel/total
known/total


known_fsm <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity > 99,
         Sqanti_structural_category %in% "full-splice_match"
  ) %>% dplyr::count()


all_fsm <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity > 0,
         Sqanti_structural_category %in% "full-splice_match"
  ) %>% dplyr::count()

known_fsm/all_fsm

known_fsm <- transdecoder_aligned %>% 
  filter(Transdecoder_ORF_type %in% "complete",
         PBid %in% non_NMD$PBid,
         p_identity > 99,
         Sqanti_structural_category %in% "full-splice_match"
  ) %>% dplyr::count()
