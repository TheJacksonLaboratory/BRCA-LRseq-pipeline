library(seqinr)
library(dplyr)

# 1. Filtering of PacBio OPEN READ FRAMES on QC passing trascripts
# 2. PARSE TRANSDECODER RESULTS to create transdecoder_map TABLE
# mapping transcripts to proteins to integrate Transdecoder, Blastp, Uniprot and Sqanti results


# Sqanti transcript annotation QC pass
sq_QCpass <- read.table("./QC_passed/Sqanti_annotation_QC_pass_transcripts.txt",
                        sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)


pep_file <- "./transdecoder/PacBio_Breast_cancer_allTranscripts_corrected.fasta.transdecoder.pep" # protein sequences
blastp_file <- "./transdecoder/blastp.out"
out_transdecoder <- "./transdecoder/Transdecoder_to_Pacbio.Rd"

pep <- read.fasta(file = pep_file, seqtype = "AA", forceDNAtolower = FALSE)
blastp_out <- read.table(blastp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)


# ****** Select peptides of transcripts passing QC
pep_names <- getName(pep)
pep_sequence <- getSequence(pep)
pep_ids <- getName(pep)
pb_ids <- gsub(".p[0-9]", "", pep_ids, perl = TRUE)
head(pep_sequence)
head(pb_ids)

# Select only peptides from QC passing transcripts *** pep_QCpass_transcripts **
pass_QC <- pb_ids %in% sq_QCpass$isoform
sum(pass_QC)
pep_QCpass_transcripts <- pep[pass_QC]
annotation <- getAnnot(pep_QCpass_transcripts)
pb_ids_pep_QCpass <- gsub(".p[0-9]", "", getName(pep_QCpass_transcripts), perl = TRUE)

# Re-write PacBio FASTA with uppercase and removing * symbol
pep_sequences <- getSequence(pep_QCpass_transcripts)
pep_modified <- lapply(seq_along(pep_sequences), function(i){
  seq_a <- paste(toupper(as.character(unlist(pep_sequences[i]))), collapse = "") 
  seq_a <- gsub("\\*", "", seq_a)
})

write.fasta(sequences = pep_modified, names = getName(pep_QCpass_transcripts), 
            file.out = "./transdecoder/QC_pass_transcripts.transdecoder.fasta.pep",
            as.string = TRUE)



has_blastp <- getName(pep_QCpass_transcripts) %in% blastp_out$V1
sum(has_blastp)

## *********** ORF type using pep_with_blastp
ORF_type <- rep("NA", length(pep_QCpass_transcripts))
fl_index <- grep("type:complete",annotation,ignore.case=T)
fivePrime_index <- grep("type:5prime_partial",annotation,ignore.case=T)
threePrime_index <- grep("type:3prime_partial",annotation,ignore.case=T)
internal_index <- grep("type:internal",annotation,ignore.case=T)
ORF_type[fl_index] <- "complete"
ORF_type[fivePrime_index] <- "5prime_partial"
ORF_type[threePrime_index] <- "3prime_partial"
ORF_type[internal_index] <- "internal"
table(ORF_type)

## *********** Blastp match
uniprot_id <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], ",", perl = TRUE)
  res2 <- strsplit(res[[1]][3], "\\|")
  uid <- res2[[1]][1]
  
})

#  Local blastp percent identity (within the region aligned in the subject)
pident <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], ",", perl = TRUE)
  res2 <- strsplit(res[[1]][3], "\\|")
  res2[[1]][2]
  
})

# Blastp E-value is the number of expected hits of similar quality (score) that could be found just by chance.
evalue <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], ",", perl = TRUE)
  res2 <- strsplit(res[[1]][3], "\\|")
  res3 <- strsplit(res2[[1]][3], " ")
  as.numeric(res3[[1]][1])
  
})

head(evalue[!has_blastp])

# Transdecoder score
tscore <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], "score=")
  res2 <- strsplit(res[[1]][2], ",")
  as.numeric(res2[[1]][1])
  
})

head(tscore[has_blastp])

# Transdecoder ORF length
tlength <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], "len:")
  res2 <- strsplit(res[[1]][2], " ")
  as.numeric(res2[[1]][1])
  
})
head(tlength[!has_blastp])

# Transdecoder ORF strand
tstrand <- lapply(annotation, function(annot){
  
  res <- strsplit(annot[[1]], ",")
  if(grepl("\\(\\+\\)", res[[1]][1], perl = TRUE)){
    return("+")
  }else{
    return("-")
  }
  
})
head(tstrand[has_blastp])

# Sqanti annotation
idx_sq <- match(pb_ids_pep_QCpass, sq_QCpass$isoform)
sum(is.na(idx_sq))

Sqanti_associated_gene <- sq_QCpass$associated_gene[idx_sq]
Sqanti_associated_gene_simp <- lapply(sq_QCpass$associated_gene[idx_sq], function(id){
  strsplit(as.character(id), "\\.")[[1]][1]
})
Sqanti_associated_gene_simp <- unlist(Sqanti_associated_gene_simp)

# ****** Annotate with Gene type with Gencode ************
gencode <- import.gff("/Users/dveiga/tools/Gencode/gencode.v30.annotation.sorted.gff.gz")
gencode <- gencode[gencode$type %in% "gene", ]

idx_gencode <- match(Sqanti_associated_gene, gencode$gene_id)

Gencode_gene_type <- gencode$gene_type[idx_gencode]
Gencode_gene_name <- gencode$gene_name[idx_gencode]
Gencode_gene_id <- gencode$gene_id[idx_gencode]
sum(is.na(Gencode_gene_id)) #novel genes, or transcripts with multiple ENSGs


# ****** Annotate with Uniprot to Ensembl Gene id ************
uniprot2ENSG <- read.table(file = "/Users/dveiga/tools/Uniprot/release_2019_04/reference_proteome/UP000005640_9606.idmapping.Ensembl.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = FALSE)

uniprot_id_modified <- lapply(unlist(uniprot_id), function(id){
  strsplit(as.character(id), "-")[[1]][1]
})
uniprot_id_modified <- unlist(uniprot_id_modified)
id_mapping <- match(uniprot_id_modified, uniprot2ENSG$V1)
Ensembl_geneID <- uniprot2ENSG$V3[id_mapping]
sum(is.na(Ensembl_geneID)) #NA correspond to protein sequences that are not assigned to gene loci
Uniprot2Ensembl_GeneID <- Ensembl_geneID



transdecoder_map <- data.frame(PBid = pb_ids_pep_QCpass,
                               PeptideId = getName(pep_QCpass_transcripts),
                               Blastp_match = unlist(uniprot_id),
                               Blastp_evalue = unlist(evalue),
                               Blastp_Uniprot2Ensembl = Ensembl_geneID,
                               Transdecoder_ORF_type= ORF_type,
                               Transdecoder_score = unlist(tscore),
                               Transdecoder_ORF_length = unlist(tlength),
                               Transdecoder_ORF_strand= unlist(tstrand),
                               Sqanti_structural_category = sq_QCpass$structural_category[idx_sq],
                               Sqanti_associated_gene_simp = Sqanti_associated_gene_simp, #Ensembl gene based on SQANTI,
                               Sqanti_isoform = sq_QCpass$isoform[idx_sq],
                               Sqanti_associated_transcript = sq_QCpass$associated_transcript[idx_sq],
                               Sqanti_chrom = sq_QCpass$chrom[idx_sq],
                               Sqanti_strand = sq_QCpass$strand[idx_sq],
                               Gencode_gene_name = Gencode_gene_name,
                               Gencode_gene_type = Gencode_gene_type,
                               stringsAsFactors = FALSE
)

head(transdecoder_map)

save(transdecoder_map, file = out_transdecoder)

# NOT applying this filter anymore
# ****** Filter #2 Select ORFs in the correct loci (Sqanti annotation) **********
# Conflict cases include NMD and retained intron transcripts that transdecoder finds a partial blastp similarity in other genes/chromosomes
# transdecoder_coding_correctLoci <- dplyr::filter(transdecoder_map,
#                                                  Sqanti_associated_gene_simp == Blastp_Uniprot2Ensembl)
# 
# transdecoder_coding_conflict <- dplyr::filter(transdecoder_map,
#                                               Sqanti_associated_gene_simp != Blastp_Uniprot2Ensembl)

nrow(transdecoder_coding_correctLoci)
nrow(transdecoder_coding_conflict)
