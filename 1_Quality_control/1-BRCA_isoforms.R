library(rtracklayer)
library(dplyr)

# QC filtering of PacBio transcripts

# Sqanti transcript annotation
sq <- read.table("./sqanti/PacBio_Breast_cancer_allTranscripts_classification.txt",
                 sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)
head(sq)

# Sqanti2 transcript annotation
# Sqanti2 has different classification from the original tool, especially in the subcategories and no fusion annotation
# Use Sqanti2 for determining junction coverage (Intropolis dataset)
# Use Sqanti2 for 3'end analysis
sq2 <- read.table("./sqanti2/PacBio_Breast_cancer_allTranscripts_classification.txt",
                 sep = "\t", stringsAsFactors = FALSE, quote = NULL, header = TRUE)


# Add min_cov from Sqanti2
idx_mapping <- match(sq$isoform, sq2$isoform)
sq$sqanti2_isoform <- sq2$isoform[idx_mapping]
all.equal(sq$isoform, sq$sqanti2_isoform)
sq$min_cov_intropolis <- sq2$min_cov[idx_mapping]

# Gffcompare annotation 
# Intron retention class codes "m", "n", "i", "y"
gffcompare_gtf <- import.gff("./gffcompare/PacBio_Breast_cancer_allTranscripts_corrected_gffcompare.annotated.gtf.gz")
gffcompare_gtf_intronRet <- gffcompare_gtf[gffcompare_gtf$class_code %in% c("m", "n", "i", "y"), ]
length(unique(gffcompare_gtf_intronRet$transcript_id))

# Add number PacBio samples where transcript was detected
# Create PacBio FL counts matrix for breast samples
# To be used in the coverage filter
sample_annot <- read.table("/Users/dveiga/analysis/Grants/Cancer-Sep2017-resubmission/GTF_stats_CancerOctober13/Samples_tissue_annotated_final.txt",
                           quote = NULL, sep = "\t", fill = TRUE, header = TRUE)

brca <- sample_annot %>% filter(Broad_classification %in% c("Breast Cancer", "Normal Breast"))

iso2flCounts <- read.table("/Users/dveiga/analysis/Grants/Cancer-Sep2017-resubmission/GTF_stats_CancerOctober13/all_samples.chained_count.txt",
                           quote = NULL, sep = "\t", header = TRUE)
rownames(iso2flCounts) <- iso2flCounts$superPBID
idx_columns <- match(brca$Sample, colnames(iso2flCounts))
idx_rows <- match(sq$isoform, iso2flCounts$superPBID)
sum(is.na(idx_rows))
sum(is.na(idx_columns))

iso2flCounts_brca <- iso2flCounts[idx_rows, idx_columns]
head(iso2flCounts_brca[, 1:57])
dim(iso2flCounts_brca)

FL_counts <- rep(NA, nrow(iso2flCounts_brca))
FL_samples <- rep(NA, nrow(iso2flCounts_brca))

for (i in 1:nrow(iso2flCounts_brca)) {
  FL_counts[i] <- sum(as.numeric(iso2flCounts_brca[i, ]), na.rm = TRUE)
  FL_samples[i] <- sum(!is.na(iso2flCounts_brca[i,]))
}

range(FL_counts)
median(FL_counts)
range(FL_samples)
median(FL_samples)

PB_coverage <- data.frame(isoform = rownames(iso2flCounts_brca),
                          FL_counts = FL_counts,
                          FL_samples = FL_samples,
                          stringsAsFactors = FALSE)
head(PB_coverage)

# Add min_cov from Sqanti2
idx_mapping2 <- match(sq$isoform, PB_coverage$isoform)
all.equal(sq$isoform, PB_coverage$isoform[idx_mapping2])
sq$FL_counts <- PB_coverage$FL_counts[idx_mapping2]
sq$FL_samples <- PB_coverage$FL_samples[idx_mapping2]


# **************** Assess 3' end reliability of FSM, ISM, NIC and NNC transcripts

# Uses SQANTI 2 ANNOTATION
# FSM and ISM categories have associated transcript id (ENST)
# Use dist_to_TTS -> distance of isoform 3' end to reference annotated end site. 
# Negative value means query ends upstream of reference.

# NIC and NNC categories have associated gene id (ENSG)
# Use diff_to_gene_TSS -> distance of isoform 3' end to the closest end of any transcripts of the matching gene
# Negative value means query ends upstream of reference.


last_pac_exons_polyA <- import.gff("./polyA_annotation/Pacbio_last_exons_polyA_support.gtf")

# Parameter sets the maximum allowed distance to an annotated 3' end for not checking intrapriming
max_dist_annotatedTTS <- 100

# Parameter determines the fraction of genomic 'A's above which the isoform will be filtered
# perc_A_downstream_TTS
a_fraction <- 80

# Rule: if dist to TTS is more than max_dist_annotatedTTS 
#       and perc_A_downstream_TTS > a_fraction, 
#       and isoform last exon does not have polyA site support
#       then => 3'end is unreliable because of possible intra-priming

# Equivalent to: Pacbio last exon does not overlap a Gencode last exon (allowing for 50bp)
#                and there is no polyA support for the last exon


# FSM category
FSM_transcripts_flagged_3end <- sq2 %>% 
  dplyr::filter(structural_category == "full-splice_match",
                !isoform %in% last_pac_exons_polyA$transcript_id,
                abs(diff_to_TTS) > max_dist_annotatedTTS,
                perc_A_downstream_TTS > a_fraction)

dim(FSM_transcripts_flagged_3end)
head(FSM_transcripts_flagged_3end)
table(FSM_transcripts_flagged_3end$subcategory)
range(abs(FSM_transcripts_flagged_3end$diff_to_TTS))
mean(abs(FSM_transcripts_flagged_3end$diff_to_TTS))

# ISM category
ISM_transcripts_flagged_3end <- sq2 %>% 
  dplyr::filter(structural_category == "incomplete-splice_match",
                !isoform %in% last_pac_exons_polyA$transcript_id,
                abs(diff_to_TTS) > max_dist_annotatedTTS,
                perc_A_downstream_TTS > a_fraction)

dim(ISM_transcripts_flagged_3end)
head(ISM_transcripts_flagged_3end)
table(ISM_transcripts_flagged_3end$subcategory)
range(abs(ISM_transcripts_flagged_3end$diff_to_TTS))
mean(abs(ISM_transcripts_flagged_3end$diff_to_TTS))

# NIC category - use diff_to_geneTTS
NIC_transcripts_flagged_3end <- sq2 %>% 
  dplyr::filter(structural_category == "novel_in_catalog",
                !isoform %in% last_pac_exons_polyA$transcript_id,
                abs(diff_to_gene_TTS) > max_dist_annotatedTTS,
                perc_A_downstream_TTS > a_fraction)

dim(NIC_transcripts_flagged_3end)
head(NIC_transcripts_flagged_3end)
table(NIC_transcripts_flagged_3end$subcategory)
range(abs(NIC_transcripts_flagged_3end$diff_to_gene_TTS))
mean(abs(NIC_transcripts_flagged_3end$diff_to_gene_TTS))

# NNC category - - use diff_to_geneTTS
NNC_transcripts_flagged_3end <- sq2 %>% 
  dplyr::filter(structural_category == "novel_not_in_catalog",
                !isoform %in% last_pac_exons_polyA$transcript_id,
                abs(diff_to_gene_TTS) > max_dist_annotatedTTS,
                perc_A_downstream_TTS > a_fraction)

dim(NNC_transcripts_flagged_3end)
head(NNC_transcripts_flagged_3end)
table(NNC_transcripts_flagged_3end$subcategory)
range(abs(NNC_transcripts_flagged_3end$diff_to_gene_TTS))
mean(abs(NNC_transcripts_flagged_3end$diff_to_gene_TTS))


# ***************** Final filtering to obtain isoforms

# FSM category
FSM_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "full-splice_match",
                exons > 1)

dim(FSM_transripts_passedQC)
nrow(FSM_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "full-splice_match" ) %>% count())

# ISM category
ISM_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "incomplete-splice_match",
                exons > 1,
                !isoform %in% ISM_transcripts_flagged_3end$isoform)

dim(ISM_transripts_passedQC)
nrow(ISM_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "incomplete-splice_match" ) %>% count())


# NIC category
cov_pass_NIC <- rbind( dplyr::filter(sq, structural_category == "novel_in_catalog", FL_samples >= 3),
                       dplyr::filter(sq, structural_category == "novel_in_catalog",min_cov_intropolis >= 5))

NIC_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "novel_in_catalog",
                exons > 1,                                           # mono-exon
                !isoform %in% NIC_transcripts_flagged_3end$isoform,  # 3'end 
                isoform %in% cov_pass_NIC$isoform,        # Coverage
                !isoform %in% gffcompare_gtf_intronRet$transcript_id) # no IR


dim(NIC_transripts_passedQC)
nrow(NIC_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "novel_in_catalog" ) %>% count())


# NNC category
cov_pass_NNC <- rbind( dplyr::filter(sq, structural_category == "novel_not_in_catalog", FL_samples >= 3),
                       dplyr::filter(sq, structural_category == "novel_not_in_catalog",min_cov_intropolis >= 5))

NNC_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "novel_not_in_catalog",
                exons > 1,
                !isoform %in% NNC_transcripts_flagged_3end$isoform,
                isoform %in% cov_pass_NNC$isoform,
                !isoform %in% gffcompare_gtf_intronRet$transcript_id,
                RTS_stage %in% "FALSE",
                all_canonical %in% "canonical")

dim(NNC_transripts_passedQC)
nrow(NNC_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "novel_not_in_catalog" ) %>% count())


# Intergenic category
cov_pass_intergenic <- rbind( dplyr::filter(sq, structural_category == "intergenic", FL_samples >= 3),
                       dplyr::filter(sq, structural_category == "intergenic",min_cov_intropolis >= 5))

intergenic_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "intergenic",
                exons > 1,
                isoform %in% cov_pass_intergenic$isoform,
                RTS_stage %in% "FALSE",
                all_canonical %in% "canonical")


dim(intergenic_transripts_passedQC)
nrow(intergenic_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "intergenic" ) %>% count())


# Antisense category
# Majority has non-canonical junctions
cov_pass_AS <- rbind( dplyr::filter(sq, structural_category == "antisense", FL_samples >= 3),
                              dplyr::filter(sq, structural_category == "antisense", min_cov_intropolis >= 5))

AS_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "antisense",
                exons > 1,
                isoform %in% cov_pass_AS$isoform,
                RTS_stage %in% "FALSE",
                all_canonical %in% "canonical")

dim(AS_transripts_passedQC)
nrow(AS_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "antisense" ) %>% count())

# Fusion category (read-through transcripts, multiple ENSG)
cov_pass_fusion <- rbind( dplyr::filter(sq, structural_category == "fusion", FL_samples >= 3),
                      dplyr::filter(sq, structural_category == "fusion", min_cov_intropolis >= 5))

Fusion_transripts_passedQC <- sq %>% 
  dplyr::filter(structural_category == "fusion",
                exons > 1,
                isoform %in% cov_pass_fusion$isoform,
                RTS_stage %in% "FALSE",
                all_canonical %in% "canonical")

dim(Fusion_transripts_passedQC)
nrow(AS_transripts_passedQC) / (
  dplyr::filter(sq, structural_category == "fusion" ) %>% count())


# Sqanti annotation for QC pass transcripts
final_set_transcripts <- rbind(FSM_transripts_passedQC,
                               ISM_transripts_passedQC,
                               NIC_transripts_passedQC,
                               NNC_transripts_passedQC,
                               intergenic_transripts_passedQC,
                               AS_transripts_passedQC,
                               Fusion_transripts_passedQC)
nrow(final_set_transcripts)

write.table(final_set_transcripts, file = "./QC_passed/Sqanti_annotation_QC_pass_transcripts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# GTF generation
# Start with Gffcompare GTF because it has exon numbers

pac_gtf <- import.gff("./gffcompare/PacBio_Breast_cancer_allTranscripts_corrected_gffcompare.annotated.gtf.gz")
pac_gtf

genc_gtf <- import.gff("~/tools/Gencode/gencode.v30.annotation.sorted.gff.gz")

# Modifying GTF columns in the PacBio GTF
# Add gene_id (Sqanti) and gene_name (pull from Gencode)
pac_gtf$transcript_name <- pac_gtf$transcript_id

idx_sq <- match(pac_gtf$transcript_id, sq$isoform)
sum(is.na(idx_sq))
pac_gtf$gene_id <- sq$associated_gene[idx_sq]
pac_gtf$structural_category <- sq$structural_category[idx_sq]
pac_gtf$associated_transcript <- sq$associated_transcript[idx_sq]
pac_gtf$gffcompare_gene_name <- pac_gtf$gene_name #this column is actually gffcompare gene name

idx_gencode <- match(pac_gtf$gene_id, genc_gtf$gene_id)
pac_gtf$gene_name <- genc_gtf$gene_name[idx_gencode]
pac_gtf$gene_type <- genc_gtf$gene_type[idx_gencode]

# Selecting annation from colnames(mcols(pac_gtf))
pac_gtf_simp <- pac_gtf[, c("transcript_id", "source", "type", "gene_id", "gene_name", "gene_type", 
                            "exon_number", "transcript_name", "structural_category", "associated_transcript")]
pac_gtf_simp

# Selecting QC passing transcripts
pac_gtf_simp_qcPass <- pac_gtf_simp[pac_gtf_simp$transcript_id %in% final_set_transcripts$isoform, ]

length(unique(pac_gtf_simp_qcPass$transcript_id))
nrow(final_set_transcripts)
pac_gtf_simp_qcPass

export.gff(pac_gtf_simp_qcPass, "./QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf")
