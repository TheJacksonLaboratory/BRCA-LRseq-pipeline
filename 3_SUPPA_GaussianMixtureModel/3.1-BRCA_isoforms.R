library(rtracklayer)

# find ballgown/ | grep gtf > list.txt
#for file in $(<list.txt); do cp "$file" samples_gtf; done


# Create input matrix for SUPPA
# TPM expression of transcripts

#local
# gff_file <- "./suppa/stringtie_suppaGTF/PacBio_NIC_NNC_ISM_Gencode30_merged.gtf"
# path_samples <- "./suppa/stringtie_suppaGTF/samples_gtf"
# out_file <- "./suppa/stringtie_suppaGTF/TPM_matrix.tsv"

#cluster
gff_file <- "PacBio_NIC_NNC_ISM_Gencode30_merged.gtf"
path_samples <- "./samples_gtf"
out_file <- "TPM_matrix.tsv"

suppa_gtf <- import.gff(gff_file)
gtf <- list.files(path = path_samples, 
                  full.names = TRUE)

samples <- list.files(path = path_samples)
samples <- gsub(".gtf", "", samples)

length(samples)


TPM_matrix <- matrix(data = NA, 
                     nrow = length(unique(suppa_gtf$transcript_id)), 
                     ncol = length(gtf))
rownames(TPM_matrix) <- unique(suppa_gtf$transcript_id)
colnames(TPM_matrix) <- samples
head(TPM_matrix[1:10, 1:10])
dim(TPM_matrix)

TPM_matrix <- TPM_matrix[ -1*which(is.na(rownames(TPM_matrix))),  ]
dim(TPM_matrix)

for (i in 1:length(gtf)) {
  
  cat("Sample: ", gtf[i], "\n")
  
  sample_gtf <- import.gff(gtf[i])
  gtf_txn <- sample_gtf[sample_gtf$type %in% "transcript", ]
  
  gtf_txn_filt <- gtf_txn[gtf_txn$transcript_id %in% rownames(TPM_matrix), ]
  
  
  idx_txn <- match(rownames(TPM_matrix), gtf_txn$transcript_id)
  
  
  if(sum(is.na(idx_txn))>0){ #at least one trasncript not quantified in this sample (NA)
    not_found <- "true"
  }else{
    not_found <- "false"
  }
    
  #head(rownames(TPM_matrix)) == head(gtf_txn$transcript_id[idx_txn])
  TPM_matrix[ , i] <- as.numeric(gtf_txn$TPM[idx_txn])
  
  cat("Sample: ", gtf[i], "Not_found=", not_found, "\n")
  
}

write.table(TPM_matrix, file = out_file, sep = "\t", quote = FALSE)

