library(dplyr)
library(parallel)
library(doParallel)
library(rtracklayer)

# Run in cluster using 1.2-BRCA_isoforms.pbs
# Remove duplicates from final set of transcripts
# GTF Generation and 
# Sqanti annotation for QC pass transcripts without duplicates

gtf_file <- "/projects/dveiga/analysis/git/BRCA_isoforms/QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf"

pac_gtf <- import.gff(gtf_file)
pac_exons <- pac_gtf[pac_gtf$type == "exon", ]
pac_exons_df <- as.data.frame(pac_exons)

txn_ids <- unique(pac_gtf$transcript_id)

res <- mclapply(1:length(txn_ids), function(i){
  
  gid <- pac_exons_df %>% filter(transcript_id == txn_ids[i]) %>% pull(gene_id) %>% unique()
  gene_txn <- pac_exons_df %>% filter(gene_id == gid) %>% pull(transcript_id) %>% unique()
  
  other_txn <- setdiff(gene_txn, txn_ids[i])
  other_txn_match <- c()
  
  if(length(other_txn) == 0){
    return(data.frame(#transcript_id = txn_ids[i],
                      txn_match = paste(c(txn_ids[i], other_txn_match), collapse = ","),
                      stringsAsFactors = FALSE)
    )
  }
  
  exons_target <- pac_exons_df %>% filter(transcript_id == txn_ids[i]) %>%
    dplyr::select(c("seqnames", "start", "end", "strand"))
  
  for (t in 1:length(other_txn)) {
    exons_subj <- pac_exons_df %>% filter(transcript_id == other_txn[t]) %>%
      dplyr::select(c("seqnames", "start", "end", "strand"))
    
    if(dplyr::all_equal(exons_target, exons_subj) == TRUE){
      other_txn_match <- c(other_txn_match, other_txn[t])  
    }
    
  }
  
  cat("Transcript ", i, "\n")
  
  return(data.frame(#transcript_id = txn_ids[i],
             txn_match = paste(sort(c(txn_ids[i], other_txn_match)), collapse = ","),
             stringsAsFactors = FALSE)
  )
  
    
}, mc.cores = 8)

df_res <- do.call(rbind, res)

# Remove duplicates
nondup <- dplyr::distinct(df_res)
nondup_ids <- rep(NA, nrow(nondup))

dim(nondup)

# Select Pbid for non-unique transcripts
for (i in 1:nrow(nondup)) {
  nondup_ids[i] <- strsplit(nondup$txn_match[i], ",")[[1]][1]  
}

pb_dup_mapping <- cbind(nondup, nondup_ids)
save(pb_dup_mapping, file = "./QC_passed/pb_ids_noDup_mapping.Rd")

load("./QC_passed/pb_ids_noDup_mapping.Rd")
dim(pb_dup_mapping)
head(pb_dup_mapping, 30)

# GTF Generation and 
# Sqanti annotation for QC pass transcripts without duplicates
pac_gtf_unique <- pac_gtf[pac_gtf$transcript_id %in% pb_dup_mapping$nondup_ids, ] 
export.gff(pac_gtf_unique, "./QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts_unique.gtf")

sq <- read.table("./QC_passed/Sqanti_annotation_QC_pass_transcripts.txt", header = TRUE, quote = NULL,sep = "\t")

sq_nonDup <- sq %>% filter(isoform %in% pb_dup_mapping$nondup_ids)
dim(sq_nonDup)

write.table(sq_nonDup, file = "./QC_passed/Sqanti_annotation_QC_pass_transcripts_unique.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
