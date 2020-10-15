library(seqinr)

# **** FASTA file for QC pass transcripts ************

cdna_seq <- read.fasta(file = "sqanti/PacBio_Breast_cancer_allTranscripts_corrected.fasta.gz",
                       seqtype = "DNA")
length(cdna_seq)
cdna_seq[1]

# Sqanti transcript annotation QC pass
sq_QCpass <- read.table("./QC_passed/Sqanti_annotation_QC_pass_transcripts_unique.txt",
                        sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)

head(sq_QCpass)

idx_keep <- getName(cdna_seq) 

cdna_seq_QCpass <- cdna_seq[getName(cdna_seq) %in% sq_QCpass$isoform]
length(cdna_seq_QCpass)

write.fasta(cdna_seq_QCpass, getName(cdna_seq_QCpass), 
            file.out = "./QC_passed/PacBio_Breast_cancer_QC_pass_transcripts_unique.fasta")

