library(dplyr)
library(seqinr)
library(doParallel)
library(parallel)

# Isoform functional annotation
# 2. Test presence/absence of transmembrane domains (TMHMM) comparing PacBio ORF and UniProt by splice variants

# INPUTS local
# input_Pacbio_to_ORF <- "./transdecoder/Transdecoder_map_aligned_annot.Rd"
# input_tmhmm_pacbio <- "./tmhmm/out_tmhmm_regions.txt"
# input_tmhmm_uniprot <- "./transdecoder/uniprot_iso_annot/out_tmhmm_regions.txt"
# 
# output_file <- "./transdecoder/uniprot_iso_annot/df_tmhmm_effect.Rd"


# INPUTS cluster
input_Pacbio_to_ORF <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment/Transdecoder_map_aligned.Rd"
input_tmhmm_pacbio <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder/tmhmm_regions.txt"
input_tmhmm_uniprot <- "/projects/dveiga/analysis/git/ts3_analysis/uniprot_iso_annot/out_tmhmm_regions.txt"

output_file <- "./df_tmhmm_effect.Rd"


# ***** PacBio ORFeome annotation
load(input_Pacbio_to_ORF)
head(transdecoder_aligned)


# ***** TMHMM regions on Pacbio proteins
pb_tmhmm <- read.table(file = input_tmhmm_pacbio, 
                       header = F, stringsAsFactors = F)

colnames(pb_tmhmm) <- c("query_name", "TMHMM_version",
                        "FT", "start", "end")
head(pb_tmhmm)

up_tmhmm <- read.table(file = input_tmhmm_uniprot, 
                       header = F, stringsAsFactors = F)

colnames(up_tmhmm) <- c("query_name", "TMHMM_version",
                        "FT", "start", "end")

uniprot_names_modified <- lapply(seq_along(up_tmhmm$query_name), function(i){
  strsplit(up_tmhmm$query_name[i], "\\|")[[1]][2]
})
uniprot_names_modified <- unlist(uniprot_names_modified)
up_tmhmm$query_name_simp <- uniprot_names_modified


res <- mclapply(seq(1:nrow(transdecoder_aligned)), function(i){
  
  hits_pb <- pb_tmhmm %>% filter(query_name %in% transdecoder_aligned$PeptideId[i], FT %in% "TMhelix") %>% nrow()
  hits_up <- up_tmhmm %>% filter(query_name_simp %in% transdecoder_aligned$Blastp_match[i], FT %in% "TMhelix") %>% nrow()
  
  if(hits_pb == 0 && hits_up == 0){
    return(data.frame()) 
  }
  
  if(hits_pb > 0 && hits_up > 0){
    tm_effect <- "No_change"
  }
  
  if(hits_pb > 0 && hits_up == 0){
    tm_effect <- "Gain"
  }
  
  if(hits_pb == 0 && hits_up > 0){
    tm_effect <- "Loss"
  }
  
  return(data.frame(pep_id = as.character(transdecoder_aligned$PeptideId[i]), 
                    effect = tm_effect,
                    splice_variant = as.character(transdecoder_aligned$Sqanti_structural_category[i]))) 
  
}, mc.cores = detectCores())


df_tmhmm_effect <- do.call(rbind, res)
dim(df_tmhmm_effect)
df_tmhmm_effect %>% group_by(splice_variant, effect) %>% dplyr::count()

save(df_tmhmm_effect, file = output_file)
