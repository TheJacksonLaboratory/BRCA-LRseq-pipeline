#!/bin/bash
#SBATCH -J transdecoder_downstream
#SBATCH -D /projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment
#SBATCH -N 1
#SBATCH -n=32
#SBATCH --mem 32G
#SBATCH -t 0-12:00:00
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -o %x.out
#SBATCH -e %x.out
set -u
dir_r=/projects/banchereau-lab/tools/3_4_4/bin
$dir_r/R CMD BATCH /projects/dveiga/analysis/git/BRCA_isoforms/transdecoder_globalAlignment/2.1-BRCA_isoforms.R

