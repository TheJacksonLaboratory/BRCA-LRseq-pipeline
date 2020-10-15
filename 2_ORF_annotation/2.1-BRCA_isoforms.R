library(seqinr)
library(Biostrings)
library(doParallel)
library(parallel)


# ***** RUN on helix
# GLOBAL PAIRWISE ALIGNMENT WITH UNIPROT
# Analysis to compute ORF Percent identity to known UNIPROT isoforms
# Align ORF from FL isoform to closest Uniprot (defined by blastp agains Uniprot)

#MODIFY PATHS
in_transdecoder <- "./Transdecoder_to_Pacbio.Rd"
in_uniprot <- "/projects/banchereau-lab/tools/Uniprot/release_2019_04/reference_proteome/UP000005640_9606_allProteins_combined.fasta"
in_pep <- "QC_pass_transcripts.transdecoder.fasta.pep"
output_aligned = "./Transdecoder_map_aligned.Rd"

# in_transdecoder <- "./transdecoder/Transdecoder_to_Pacbio.Rd"
# in_uniprot <- "/Users/dveiga/tools/Uniprot/release_2019_04/reference_proteome/UP000005640_9606_allProteins_combined.fasta"
# in_pep <- "./transdecoder/QC_pass_transcripts.transdecoder.fasta.pep"
# output_aligned = "./transdecoder/Transdecoder_map_aligned.Rd"


load(file = in_transdecoder) #load transdecoder_map
head(transdecoder_map)
uniprot_fa <- read.fasta(file = in_uniprot)
length(uniprot_fa)
uniprot_names <- getName(uniprot_fa)
uniprot_sequences <- getSequence(uniprot_fa)

# Transdecoder ORF sequences
pep <- read.fasta(file = in_pep)
length(pep)
annotation <- getAnnot(pep)
pep_names <- getName(pep)
pep_sequences <- getSequence(pep)


pid_align <- mclapply(seq_along(transdecoder_map$PBid), function(i){
  
  idx_a <- match(transdecoder_map$PeptideId[i], pep_names )
  idx_b <- grep(transdecoder_map$Blastp_match[i], uniprot_names)
  
  seq_a <- paste(toupper(as.character(unlist(pep_sequences[idx_a]))), collapse = "") 
  #seq_a == transdecoder_map$AAseq[i]
  seq_a <- gsub("\\*", "", seq_a)
  seq_b <- paste(toupper(as.character(unlist(uniprot_sequences[idx_b]))), collapse = "")
  
  if(nchar(seq_a) == 0){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -1)) 
    #code for query sequence is length 0
  }
  
  if(nchar(seq_b) == 0){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -2)) 
    #code for subject sequence is length 0
  }
  
  if(any(!strsplit(seq_a, "")[[1]] %in% names(AMINO_ACID_CODE))){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -3)) 
    #code for query/subject sequence contain non-allowed characters
  }
  
  if(any(!strsplit(seq_b, "")[[1]] %in% names(AMINO_ACID_CODE))){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -3)) 
    #code for query/subject sequence contain non-allowed characters
  }
  
  # # pairwiseAlignment() cannot be used on very long sequences (~ 50kb)
  # if( (nchar(seq_a) * nchar(seq_b)) > .Machine$integer.max ){
  #   return(-4) #code for aligment is too long
  # }
  
  #Pairwise global alignment of protein sequences using the Needleman-Wunsch algorithm
  an.error.occured <- FALSE
  tryCatch( {palign <- pairwiseAlignment(AAString(seq_a), AAString(seq_b),
                                         substitutionMatrix = "BLOSUM50" ,
                                         gapOpening = 10, gapExtension = 4,
                                         type = "global") },
            error = function(e) {an.error.occured <<- TRUE} )
  #print(an.error.occured)
  if(an.error.occured){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -5)) 
    #code for pairwiseAlignment() fail
  }
  
  #Percent identity of alignments ("PID1" account for aligned positions and gaps)
  an.error.occured <- FALSE
  tryCatch( {pid_result <- pid(palign, type = "PID1") },
            error = function(e) {an.error.occured <<- TRUE} )
  
  if(an.error.occured){
    return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = -6)) 
    #code for pid() fail
  }
  
  return(list(pep_id = as.character(transdecoder_map$PeptideId[i]), value = pid_result)) 
  
}, mc.cores = 32)


df <- data.frame(matrix(unlist(pid_align), nrow=length(pid_align), byrow=T), stringsAsFactors = F)
colnames(df) <- names(pid_align[[1]])
df$value <- round(as.numeric(df$value), 2)
head(df)

idx_pid <- match(transdecoder_map$PeptideId, df$pep_id)
sum(as.character(transdecoder_map$PeptideId) != df$pep_id[idx_pid])
transdecoder_map$p_identity <- df$value[idx_pid]
transdecoder_aligned <- transdecoder_map

save(transdecoder_aligned, file = output_aligned)

# code -5 => code for pairwiseAlignment() fail, unknown reason few cases
# code -2 => no blastp match
table(transdecoder_aligned$p_identity[which(transdecoder_aligned$p_identity < 0)])

