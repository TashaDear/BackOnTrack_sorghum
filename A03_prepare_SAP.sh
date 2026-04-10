#!/bin/bash
#SBATCH --partition short
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --account=backontrack
#SBATCH -t 12:0:0
#SBATCH --mem=200g
#SBATCH --job-name WP3_prepare_SAP
#############################################################
source /home/tasha/miniconda3/etc/profile.d/conda.sh

#Paths
GFF_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"
FA_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa"

#Input 
PATH_SAP="/home/tasha/backontrack/input/SAP/"

#Output
PATH_SAP_out="/home/tasha/backontrack/input/vcf_annotations/SAP/"

#############################################################
#A) Annotate degeneracy
#############################################################

conda activate R_env

R CMD BATCH --max-ppsize=500000 Annotate_SAP.R 

#############################################################
#B) Prepare file
#############################################################

conda activate vcf_bcf

#1) Filter out multi-allelic sites
bcftools view -m2 -M2 \
"${PATH_SAP_out}SAP.imputed_upd_CDS_deg.vcf.gz" \
-Ov -o "${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic.vcf.gz"

#index
tabix -p vcf -f "${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic.vcf.gz"

#2) Annotate with bcftools
bcftools csq -f "${FA_FILE}" \
-g "${GFF_FILE}" \
"${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic.vcf.gz" \
-o "${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic_annotated.vcf.gz" \
-O z

#index
tabix -p vcf -f "${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic_annotated.vcf.gz" 

#3) Extract the necessary information
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Degeneracy\t%INFO/BCSQ\n' \
  "${PATH_SAP_out}SAP.imputed_upd_CDS_deg_biallelic_annotated.vcf.gz" | \
awk -F'\t' 'BEGIN {OFS="\t"} {
  split($6, bcsq, "|");
  print $1, $2, $3, $4, $5,
        bcsq[1], bcsq[2], bcsq[3], bcsq[4], bcsq[5], bcsq[6]
}' > "${PATH_SAP_out}SAP_CDS_biallelic_final.txt"

#############################################################
