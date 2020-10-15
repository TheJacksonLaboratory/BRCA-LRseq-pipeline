library(dplyr)
library(seqinr)
library(doParallel)
library(parallel)

# Isoform functional annotation
# 3. Test change of localization (DeepLoc) comparing PacBio ORF and UniProt by splice variants
# ***** RUN on cluster using 2.4-BRCA_isoforms.pbs

# INPUTS local
# input_Pacbio_to_ORF <- "./transdecoder/Transdecoder_map_aligned_annot.Rd"
# input_deeploc_pacbio <- "./deeploc_transdecoder/output/"
# input_deeploc_uniprot <- "./deeploc_uniprot/output/"
# output_file <- "./transdecoder/uniprot_iso_annot/df_deeploc_effect.Rd"

# INPUTS cluster
input_Pacbio_to_ORF <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment/Transdecoder_map_aligned.Rd"
input_deeploc_pacbio <- "/projects/dveiga/analysis/git/BRCA_isoforms/deeploc/pbs"
input_deeploc_uniprot <- "/projects/dveiga/analysis/git/ts3_analysis/uniprot_deeploc/pbs"
output_file <- "./df_deeploc_effect.Rd"


# ***** PacBio ORFeome annotation
load(input_Pacbio_to_ORF)
head(transdecoder_aligned)

# Read DeepLoc for PacBio
files <- list.files(path = input_deeploc_pacbio, pattern = "*.txt",
                    full.names = TRUE)
pred <- lapply(seq_along(files), function(i){
  read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
pred_pacbio <- do.call(rbind, pred)
dim(pred_pacbio)

# Read DeepLoc for UniProt
files <- list.files(path = input_deeploc_uniprot, pattern = "*.txt",
                    full.names = TRUE)
pred <- lapply(seq_along(files), function(i){
  read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = NULL)
})
pred_uniprot <- do.call(rbind, pred)
dim(pred_uniprot)

uniprot_id <- lapply(seq_along(pred_uniprot$ID), function(i){
  
  res <- strsplit(pred_uniprot$ID[i], "\\|", perl = TRUE)
  uid <- res[[1]][2]
  
})
pred_uniprot$Uniprot_simp <- unlist(uniprot_id)


res <- mclapply(seq(1:nrow(transdecoder_aligned)), function(i){
  
  hits_pb <- pred_pacbio %>% filter(ID %in% transdecoder_aligned$PeptideId[i])
  
  hits_up <- pred_uniprot %>% filter(Uniprot_simp %in% transdecoder_aligned$Blastp_match[i])
  
  #no prediction
  if(nrow(hits_pb) == 0 || nrow(hits_up) == 0){
    return(data.frame()) 
  }
  
  if(hits_pb$Location[1] == hits_up$Location[1]){
    effect <- "No_change"
  }
  
  if(hits_pb$Location[1] != hits_up$Location[1]){
    effect <- "Change"
  }
  
  
  return(data.frame(pep_id = as.character(transdecoder_aligned$PeptideId[i]), 
                    PacBio_location = hits_pb$Location[1],
                    UniProt_location = hits_up$Location[1],
                    effect = effect,
                    splice_variant = as.character(transdecoder_aligned$Sqanti_structural_category[i]))) 
  
}, mc.cores = 32)

df_deeploc_effect <- do.call(rbind, res)
dim(df_deeploc_effect)
df_deeploc_effect %>% group_by(splice_variant, effect) %>% dplyr::count()

save(df_deeploc_effect, file = output_file)
