library(rhmmer)
library(dplyr)
library(seqinr)
library(doParallel)
library(parallel)

# Isoform functional annotation
# 1. Test presence/absence of PFAM domains comparing PacBio ORF and UniProt by splice variants

# INPUTS local
# input_Pacbio_to_ORF <- "./transdecoder/Transdecoder_map_aligned_annot.Rd"
# input_pfam_pacbio <- "./transdecoder/pfam.domtblout"
# input_pfam_uniprot <- "./transdecoder/uniprot_iso_annot/pfam.domtblout"
# 
# output_file <- "./transdecoder/uniprot_iso_annot/df_pfam_effect.Rd"

# INPUTS cluster
input_Pacbio_to_ORF <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment/Transdecoder_map_aligned.Rd"
input_pfam_pacbio <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder/pfam.domtblout"
input_pfam_uniprot <- "/projects/dveiga/analysis/git/ts3_analysis/uniprot_iso_annot/pfam.domtblout"

output_file <- "./df_pfam_effect.Rd"

# ***** PacBio ORFeome annotation
load(input_Pacbio_to_ORF)
head(transdecoder_aligned)


# ***** PFAM domains on Pacbio proteins
# ievalue - a stringent measure of how reliable the domain may be
# E-10 cutoff is optimal (achieves a almost 0 change in the FSM variants)

pfam <- read_domtblout(input_pfam_pacbio)
pfam <- pfam %>% as.data.frame()
pfam_sign <- pfam %>% filter(domain_ievalue < 10E-10) 
#pfam_sign <- pfam %>% filter(domain_cevalue < 10E-5)
rm(pfam)
head(pfam_sign, 10)
length(unique(pfam_sign$domain_accession))

uprot_pfam <- read_domtblout(input_pfam_uniprot)
uprot_pfam <- uprot_pfam %>% as.data.frame()
uprot_pfam_sign <- uprot_pfam %>% filter(domain_ievalue < 10E-10)
#uprot_pfam_sign <- uprot_pfam %>% filter(domain_cevalue < 10E-5)
rm(uprot_pfam)
head(uprot_pfam_sign, 10)
length(unique(uprot_pfam_sign$domain_accession))

uniprot_names_modified <- lapply(seq_along(uprot_pfam_sign$query_name), function(i){
  strsplit(uprot_pfam_sign$query_name[i], "\\|")[[1]][2]
})
uniprot_names_modified <- unlist(uniprot_names_modified)
uprot_pfam_sign$query_name_simp <- uniprot_names_modified

res <- mclapply(seq(1:nrow(transdecoder_aligned)), function(i){
  
  dom_pb <- pfam_sign %>% filter(query_name %in% transdecoder_aligned$PeptideId[i]) %>% mutate(source = "Pacbio")
  dom_up <- uprot_pfam_sign %>% filter(query_name_simp %in% transdecoder_aligned$Blastp_match[i]) %>% mutate(source = "Uniprot")
  
  comb <- bind_rows(dom_pb, dom_up) %>% group_by(domain_accession, source) %>% dplyr::count()
  
  if(nrow(comb) == 0){
    return(data.frame()) 
  }
  
  dom_list <- unique(comb$domain_accession)
  dom_effect <- rep("NA", length(dom_list))
  pep_list <- rep(as.character(transdecoder_aligned$PeptideId[i]), length(dom_list))
  sqanti_annot <- rep(as.character(transdecoder_aligned$Sqanti_structural_category[i]), length(dom_list))
  
  #for each domain accession
  for (j in 1:length(dom_list)) {
    n_pb <- dom_pb %>% filter(domain_accession %in% dom_list[j]) %>% nrow()
    n_up <- dom_up %>% filter(domain_accession %in% dom_list[j]) %>% nrow()
    
    
    if(n_pb > 0 && n_up > 0){
      dom_effect[j] <- "No_change"
    }
    
    if(n_pb > 0 && n_up == 0){
      dom_effect[j] <- "Gain"
    }
    
    if(n_pb == 0 && n_up > 0){
      dom_effect[j] <- "Loss"
    }
    
    # if(n_pb > n_up){
    #   dom_effect[j] <- "Gain"
    # }
    # if(n_pb < n_up){
    #   dom_effect[j] <- "Loss"
    # }
    # if(n_pb == n_up){
    #   dom_effect[j] <- "No_change"
    # }
    
  }
  
  return(data.frame(pep_id = pep_list, pfam_accession = dom_list, effect = dom_effect,
                    splice_variant = sqanti_annot)) 
  
}, mc.cores = 32)


df_pfam_effect <- do.call(rbind, res)
dim(df_pfam_effect)
df_pfam_effect %>% group_by(splice_variant, effect) %>% dplyr::count()

save(df_pfam_effect, file = output_file)
