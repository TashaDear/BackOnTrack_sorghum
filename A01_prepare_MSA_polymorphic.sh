#!/bin/bash
#curl -O https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-55/variation/vep/sorghum_bicolor_vep_55_Sorghum_bicolor_NCBIv3.tar.gz
#############################################################
source /home/tasha/miniconda3/etc/profile.d/conda.sh

#Paths
GFF_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"
FA_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa"

#Input 
PATH_MSA="/home/tasha/backontrack/input/MSA/"
PATH_SAP="/home/tasha/backontrack/input/SAP/"

#Output
PATH_MSA_out="/home/tasha/backontrack/input/vcf_annotations/MSA/"
PATH_SAP_out="/home/tasha/backontrack/input/vcf_annotations/SAP/"
PATH_out="/home/tasha/backontrack/input/vcf_annotations/outputs/"

OUTGROUPS="Zea_mays_1 Erianthus_rufipilus" 

#############################################################
#A) Prepare SAP file
#############################################################

conda activate vcf_bcf

#1) Update chromosome field
bcftools annotate \
  --rename-chrs <(bcftools view -h "${PATH_SAP}SAP.imputed.vcf.gz" \
                  | grep '^##contig' \
                  | sed -E 's/##contig=<ID=([^,>]+).*/\1/' \
                  | awk '{new=$1; sub(/^Chr0*/, "", new); print $1 "\t" new}') \
  "${PATH_SAP}SAP.imputed.vcf.gz" \
  -Oz -o "${PATH_SAP_out}SAP.imputed_upd.vcf.gz"
  
#index
tabix -p vcf -f "${PATH_SAP_out}SAP.imputed_upd.vcf.gz"

#2) Extract chromosome and start/end positions for exons and CDS
gzip -dc "${GFF_FILE}" | awk -F '\t' '$3 == "CDS"  {print $1 "\t" $4 "\t" $5}' > "${PATH_SAP_out}CDS_positions.txt"

#3) Filter SAP -> keep sites in CDS
bcftools view -R "${PATH_SAP_out}CDS_positions.txt" \
"${PATH_SAP_out}SAP.imputed_upd.vcf.gz" \
-Oz -o "${PATH_SAP_out}SAP.imputed_upd_CDS.vcf.gz"

#index
tabix -p vcf -f "${PATH_SAP_out}SAP.imputed_upd_CDS.vcf.gz"

#############################################################
#B) Prepare MSA file
#############################################################

#1) Filter MSA -> keep sites in CDS
bcftools view \
-R "${PATH_SAP_out}CDS_positions.txt" \
"${PATH_MSA}full_msa.vcf.gz" \
-Oz -o "${PATH_MSA_out}MSA_CDS.vcf.gz" 

bcftools index "${PATH_MSA_out}MSA_CDS.vcf.gz"

#2) Remove outgroups not needed for the analysis
bcftools view \
-s ^Coix_aquatica,Miscanthus_sinensis_1,Miscanthus_sinensis_2,Eremochloa_ophiuroides,Zea_mays_2,Panicum_hallii,Paspalum_vaginatum \
-i 'GT[0]!="." && GT[1]!="."' \
"${PATH_MSA_out}MSA_CDS.vcf.gz" \
-Oz -o "${PATH_MSA_out}MSA_CDS_outgroups.vcf.gz" 

#3) Trim alleles  
bcftools view \
--threads 20 \
--trim-alt-alleles \
-Ou "${PATH_MSA_out}MSA_CDS_outgroups.vcf.gz" | \
bcftools annotate \
-x ^FORMAT/GT \
-Oz -o "${PATH_MSA_out}MSA_CDS_outgroups_trimmed.vcf.gz"

#4) Filter out multiallelic sites and indels in ALT or REF
bcftools view \
  --threads 20 \
  --max-alleles 2 \
  --exclude-types indels \
  --exclude 'ALT="*" || strlen(REF) > 1' \
  -Ou "${PATH_MSA_out}MSA_CDS_outgroups_trimmed.vcf.gz" | \
bcftools annotate \
  -x 'FORMAT,^FORMAT/GT' \
  -Oz -o "${PATH_MSA_out}MSA_CDS_outgroups_trimmed_no_indels.vcf.gz"
  
#index
tabix -p vcf -f "${PATH_MSA_out}MSA_CDS_outgroups_trimmed_no_indels.vcf.gz"

#5) Find intersection of sites between MSA and SAP
bcftools isec -O z -p "${PATH_MSA_out}intersected_temp_part1" \
"${PATH_MSA_out}MSA_CDS_outgroups_trimmed_no_indels.vcf.gz" \
"${PATH_SAP_out}SAP.imputed_upd_CDS.vcf.gz"

#Rename vcf file, 0003 -> records in SAP, sites shared with MSA 
cp "${PATH_MSA_out}intersected_temp_part1/0003.vcf.gz" "${PATH_SAP_out}SAP.imputed_upd_CDS_filtered.vcf.gz"

#index
tabix -p vcf -f "${PATH_SAP_out}SAP.imputed_upd_CDS_filtered.vcf.gz"

#############################################################
#C) Merge MSA with SAP
#############################################################

#1) Merge with SAP
bcftools merge \
--threads 20 \
--missing-to-ref \
"${PATH_MSA_out}MSA_CDS_outgroups_trimmed_no_indels.vcf.gz" \
"${PATH_SAP_out}SAP.imputed_upd_CDS_filtered.vcf.gz" \
-O z -o "${PATH_MSA_out}MSA_merged_SAP.vcf.gz"

#2) Filter out all sites that are multiallelic - no longer trim here 
bcftools view \
--threads 20 \
--max-alleles 2 \
-Oz \
-o "${PATH_MSA_out}MSA_merged_SAP_biallelic.vcf.gz" \
"${PATH_MSA_out}MSA_merged_SAP.vcf.gz"

# Index
tabix -p vcf -f "${PATH_MSA_out}MSA_merged_SAP_biallelic.vcf.gz"

#############################################################
#D) Annotate merged vcf
#############################################################

conda activate R_env

R CMD BATCH --max-ppsize=500000 Annotate_MSA.R

#############################################################
#E) Prepare file with polymorphic sites
#############################################################

conda activate vcf_bcf

#1) Index file
tabix -p vcf -f "${PATH_MSA_out}MSA_merged_SAP_biallelic_deg_AA.vcf.gz"

#2) Identify polymorphic sites
bcftools isec -O z -p "${PATH_MSA_out}intersected_temp_part2" \
"${PATH_MSA_out}MSA_merged_SAP_biallelic_deg_AA.vcf.gz" \
"${PATH_SAP_out}SAP.imputed_upd_CDS_filtered.vcf.gz"

#rename file, records from Ingroups, sites shared by SAP (filtered)
cp "${PATH_MSA_out}intersected_temp_part2/0002.vcf.gz" "${PATH_out}Polymorphic_sites.vcf.gz"

#3) Filter out outgroups
bcftools view \
-s ^Zea_mays_1,Erianthus_rufipilus \
"${PATH_out}Polymorphic_sites.vcf.gz" \
-Oz -o "${PATH_out}Polymorphic_sites_filtered.vcf.gz"

#index
tabix -p vcf -f "${PATH_out}Polymorphic_sites_filtered.vcf.gz"

#4) Annotate file
bcftools csq -f "${FA_FILE}" \
-g "${GFF_FILE}" \
"${PATH_out}Polymorphic_sites_filtered.vcf.gz" \
-o "${PATH_out}Polymorphic_sites_filtered_annotated.vcf.gz" \
-O z

#index
tabix -p vcf -f "${PATH_out}Polymorphic_sites_filtered_annotated.vcf.gz"

#5) Create txt file (annotations for ingroups, biallelic, CDS)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Degeneracy\t%INFO/AA\t%INFO/AA_info\t%INFO/BCSQ\n' \
"${PATH_out}Polymorphic_sites_filtered_annotated.vcf.gz" | \
awk -F'\t' 'BEGIN{OFS="\t"} {
  pmaj="NA";
  if($7 ~ /p_major_ancestral/) {
    match($7, /'\''p_major_ancestral'\''[ ]*:[ ]*([0-9.eE+-]+)/, arr);
    if(arr[1] != "") pmaj=arr[1];
  }
  split($8, bcsq, "|");
  print $1,$2,$3,$4,$5,$6,pmaj,
        (length(bcsq) > 0 ? bcsq[1] : ""),
        (length(bcsq) > 1 ? bcsq[2] : ""),
        (length(bcsq) > 2 ? bcsq[3] : ""),
        (length(bcsq) > 3 ? bcsq[4] : ""),
        (length(bcsq) > 4 ? bcsq[5] : ""),
        (length(bcsq) > 5 ? bcsq[6] : "")
}' > "${PATH_out}Polymorphic_sites_filtered_annotated_final.txt"

#############################################################
#F) Prepare data for sharing
#############################################################

#1) Intersect
bcftools isec \
  -n=2 \
  -w1 \
  -Oz \
  -o "${PATH_MSA_out}MSA_CDS_polymorphic_sites.vcf.gz" \
  "${PATH_MSA_out}MSA_CDS_outgroups_trimmed_no_indels.vcf.gz" \
  "${PATH_SAP_out}SAP.imputed_upd_CDS.vcf.gz"
  
#2) Prepare output file
bcftools query -l "${PATH_MSA_out}MSA_CDS_polymorphic_sites.vcf.gz" | \
paste -sd'\t' - | \
sed 's/^/CHROM\tPOS\tREF\tALT\t/' \
> "${PATH_out}MSA_CDS_polymorphic_sites_final.txt"

#Append
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
"${PATH_MSA_out}MSA_CDS_polymorphic_sites.vcf.gz" \
>> "${PATH_out}MSA_CDS_polymorphic_sites_final.txt"

#############################################################
#G) Sanity check
#############################################################

###Input
bcftools query -f'%CHROM:%POS %REF %ALT\n' "${PATH_SAP_out}SAP.imputed_upd_CDS_filtered.vcf.gz" > "${PATH_out}input_alleles_final.txt"

###Output
bcftools query -f'%CHROM:%POS %REF %ALT\n' "${PATH_out}Polymorphic_sites_filtered_annotated.vcf.gz" > "${PATH_out}output_alleles_final.txt"

#############################################################
