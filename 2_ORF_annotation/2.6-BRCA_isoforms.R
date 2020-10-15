library(seqinr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(doParallel)
library(parallel)


# *********** PTC Premature termination codon analysis
#

# Local variables
input_Pacbio_to_ORF <- "./transdecoder_globalAlignment/Transdecoder_map_aligned.Rd"
input_cds_transdecoder <- "./transdecoder/PacBio_Breast_cancer_allTranscripts_corrected.fasta.transdecoder.gff3"
input_Pacbio_gtf <- "QC_passed/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf"
output_file <- "./ORF_effect/df_PTC_effect.Rd"

#helix variables
input_Pacbio_to_ORF <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment/Transdecoder_map_aligned.Rd"
input_cds_transdecoder <- "/projects/dveiga/analysis/git/BRCA_isoforms/transdecoder/PacBio_Breast_cancer_allTranscripts_corrected.fasta.transdecoder.gff3"
input_Pacbio_gtf <- "/projects/dveiga/analysis/git/BRCA_isoforms/QC_passed/GTF_unzipped/PacBio_Breast_Cancer_all_QC_pass_transcripts.gtf"
output_file <- "./df_PTC_effect.Rd"



# ***** PacBio ORFeome annotation
load(input_Pacbio_to_ORF)
head(transdecoder_aligned)

# ***** PacBio CDS coordinates
cds_gff <- import.gff3(input_cds_transdecoder)
cds_gff <- cds_gff[cds_gff$type == "CDS"]
head(cds_gff)


# ORF_ids <- 
# transdecoder_aligned$PBid %in% as.vector(seqnames(cds_gff))  
# 
# ORF_ids <-
# transdecoder_aligned$PeptideId %in% as.vector(cds_gff$Parent)  
# sum(ORF_ids)

# ***** PacBio Exon sequences
pac_gtf <- import.gff(input_Pacbio_gtf)
pac_gtf_exons <- pac_gtf[pac_gtf$type == "exon", ]


res <- mclapply(seq(1:nrow(transdecoder_aligned)), function(i){

#res <- lapply(seq(1:200), function(i){  
 
  #Obtain CDS coordinates
  id_cds <- match(transdecoder_aligned$PeptideId[i], as.vector(cds_gff$Parent))
  cds_coord <- cds_gff[id_cds, ]
  
  pb_id_exons <- pac_gtf_exons[pac_gtf_exons$transcript_id == transdecoder_aligned$PBid[i]]
  exon_lengths <- width(pb_id_exons)
  
  if(as.character(strand(pb_id_exons[1])) == "-"){
    exon_lengths <- rev(exon_lengths)
  }
    
  exon_nstarts <- rep(NA, length(pb_id_exons))
  exon_nends <- rep(NA, length(pb_id_exons))
  
  # Map to novel 1-based coordinates
  exon_nstarts[1] <- 1  
  exon_nends[1] <- exon_lengths[1]
  
  for (j in 2:length(exon_lengths)) {
    exon_nstarts[j] <- exon_nends[j-1] + 1
    exon_nends[j] <- exon_nstarts[j] + exon_lengths[j] - 1
  }
  
  exon_1based <- IRanges(start = exon_nstarts,
                         end = exon_nends)
  width(exon_1based) == width(pb_id_exons)
  
  hits <- findOverlaps(IRanges(end(cds_coord), end(cds_coord)), exon_1based)
  stopCodonExon <- subjectHits(hits)
  lastExon <- length(pb_id_exons)
  
  #PTC locates > 50-55nt upstream the 3' most junction (end of before last exon)
  #Some genes have single 3'-untranslated exons and the normal stop codons in the -1 exon, but less
  # than 50bp away from the end of -1 exon
  #distance to last exon junction
  deltaStopLastEJ <- abs(end(cds_coord) - end(exon_1based[lastExon-1]))
  
  # distance to nearest downstream exon junction (splice site)
  deltaStopEJ <- abs(end(cds_coord) - end(exon_1based[stopCodonExon]))
  
  hits <- findOverlaps(IRanges(start(cds_coord), start(cds_coord)), exon_1based)
  startCodonExon <- subjectHits(hits)
  
  #PTC decision rule
  if( (stopCodonExon < lastExon) & deltaStopLastEJ > 55){
    PTC <- "yes"
  }else{
    PTC <- "no"
  }

  # 3prime_partial and internal ORF: no true stop codon found by Transdecoder
   if(transdecoder_aligned$Transdecoder_ORF_type[i] %in% c("3prime_partial", "internal")){
     PTC <- "NA"
   }
  
  
  
  return(data.frame(PBid = transdecoder_aligned$PBid[i],
                    PTC_status = PTC
          )
  ) 
  

}, mc.cores = 64) #mclapply
#}) #lapply

warnings()

df_PTC_effect <- do.call(rbind, res)
dim(df_PTC_effect)

save(df_PTC_effect, file = output_file)

# df_PTC_effect %>% filter(ORF_type %in% "complete") %>%
#   group_by(splice_variant, PTC_status) %>% dplyr::count()



