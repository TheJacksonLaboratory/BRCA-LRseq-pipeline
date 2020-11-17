library(rtracklayer)


# SUPPA GTF - 
# Create a GTF with NIC, NNC and ISM transcripts merged with Gencode  


# Filter QC passing non-redundant transcripts
sq_QCpass <- read.table("/Users/dveiga/analysis/github/BRCA_isoforms/QC_passed/Sqanti_annotation_QC_pass_transcripts_unique.txt",
                        sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)

pac_gtf <- import.gff("/Users/dveiga/analysis/github/BRCA_isoforms/QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts_unique.gff.gz")
genc_gtf <- import.gff("~/tools/Gencode/gencode.v30.annotation.sorted.gff.gz")

table(sq_QCpass$structural_category)

classes_gtf <- c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")

pb_ids <- sq_QCpass %>% filter(structural_category %in% classes_gtf) %>%
  pull(isoform)

length(pb_ids)

pac_gtf_subset <- pac_gtf[pac_gtf$transcript_id %in% pb_ids, ]
length(unique(pac_gtf_subset$transcript_id))

#Selecting annation from colnames(mcols(te_gtf))
pac_gtf_subset_filt <- pac_gtf_subset[, c("transcript_id", "source", "type", "gene_id", "exon_number", 
                                          "transcript_name", "gene_name")]

genc_gtf_filt <- genc_gtf[, c("transcript_id", "source", "type", "gene_id", "exon_number", 
                              "transcript_name", "gene_name")]

# Merge GTFs
comb_gtf <- c(pac_gtf_subset_filt, genc_gtf_filt)

length(unique(comb_gtf$transcript_id))
length(unique(pac_gtf_subset_filt$transcript_id)) + length(unique(genc_gtf_filt$transcript_id))       

head(comb_gtf)

export.gff(comb_gtf, "./suppa/stringtie_suppaGTF/PacBio_NIC_NNC_ISM_Gencode30_merged.gtf")

