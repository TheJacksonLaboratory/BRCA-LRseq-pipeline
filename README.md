# BRCA-LRseq-pipeline
This is the repository of the paper **A comprehensive long-read isoform analysis platform and sequencing resource for breast cancer**, published in [*Science Advances (2022)*](https://www.science.org/doi/10.1126/sciadv.abg6711).

The repository contains the pipeline used to analyze RNA isoforms and ORF products sequenced from breast cancer samples using PacBio long-read sequencing.


## Bulk download of long-read isoforms

- [GTF file](https://github.com/TheJacksonLaboratory/BRCA-LRseq-pipeline/blob/main/data/QC_pass/PacBio_Breast_Cancer_all_QC_pass_transcripts_unique.gff.gz) containing isoforms passing quality control
- [FASTA file](https://github.com/TheJacksonLaboratory/BRCA-LRseq-pipeline/blob/main/data/QC_pass/PacBio_Breast_cancer_QC_pass_transcripts_unique.fasta.gz) with isoform sequences
- [SQANTI annotation](https://github.com/TheJacksonLaboratory/BRCA-LRseq-pipeline/blob/main/data/QC_pass/Sqanti_annotation_QC_pass_transcripts_unique.txt.gz) of isoforms
- [Open read frames](https://github.com/TheJacksonLaboratory/BRCA-LRseq-pipeline/blob/main/data/QC_pass/QC_pass_transcripts.transdecoder.fasta.pep.zip) predicted by Transdecoder

The data and source code (pipeline V.1) can also be downloaded from Zenodo using this [link](https://doi.org/10.5281/zenodo.5449836). 
