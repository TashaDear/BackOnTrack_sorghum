#!/bin/bash
#SBATCH --partition short
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --account=backontrack
#SBATCH -t 4-00:00:00
#SBATCH --mem=200g
#SBATCH --job-name WP3_prepare_MSA 
#############################################################
source /home/tasha/miniconda3/etc/profile.d/conda.sh

#Paths
GFF_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"
FA_FILE="/home/tasha/sorghum_vep/data/references/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa"

#Input (raw)
PATH_MSA="/home/tasha/backontrack/input/MSA/"
PATH_SAP="/home/tasha/backontrack/input/SAP/"

#Output
PATH_MSA_out="/home/tasha/backontrack/input/vcf_annotations/MSA/"
PATH_SAP_out="/home/tasha/backontrack/input/vcf_annotations/SAP/"
PATH_out="/home/tasha/backontrack/input/vcf_annotations/MSA_SAP_merged/"

#############################################################
#A) Prepare file
#############################################################
conda activate vcf_bcf

#1) Copy and rename vcf file, 0000 -> records private to first file (MSA)
cp "${PATH_MSA}intersected_temp_part1/0000.vcf.gz" "${PATH_MSA}MSA_monomorphic.vcf.gz" 

#index
tabix -p vcf -f "${PATH_MSA}MSA_monomorphic.vcf.gz"

#############################################################
#B) Modify the output file
#############################################################

input_vcf="${PATH_MSA_out}MSA_monomorphic.vcf.gz"
output_vcf="${PATH_MSA_out}MSA_monomorphic_modified.vcf.gz"

zcat "$input_vcf" | awk -F'\t' 'BEGIN {OFS="\t"} {
    if (substr($0, 1, 1) == "#") {
        print
    } else {
        if ($5 == ".") {
            # Randomly generate an ALT allele if needed
            rand_alt = int(rand() * 4);
            if (rand_alt == 0) new_alt = "A";
            else if (rand_alt == 1) new_alt = "G";
            else if (rand_alt == 2) new_alt = "C";
            else if (rand_alt == 3) new_alt = "T";

            # Ensure new ALT allele is different from REF allele
            while (new_alt == $4) {
                rand_alt = int(rand() * 4);
                if (rand_alt == 0) new_alt = "A";
                else if (rand_alt == 1) new_alt = "G";
                else if (rand_alt == 2) new_alt = "C";
                else if (rand_alt == 3) new_alt = "T";
            }

            $5 = new_alt;
        }
        print
    }
}' | bgzip > "$output_vcf"

#index
tabix -p vcf -f "${PATH_MSA_out}MSA_monomorphic_modified.vcf.gz"

#############################################################
#C) Annotate with VEP
#############################################################

conda activate vep_env

#1) Annotate
vep \
  --format vcf \
  --input_file "${PATH_MSA_out}MSA_monomorphic_modified.vcf.gz" \
  --output_file "${PATH_MSA_out}MSA_monomorphic_modified_vep.txt" \
  --cache \
  --cache_version 55 \
  --dir_cache "${PATH_CACHE_GR}" \
  --offline \
  --tab \
  --force_overwrite \
  --genomes \
  --species sorghum_bicolor \
  --fasta ${FA_FILE} \
  --symbol \
  --gene_phenotype \
  --regulatory \
  --sift s \
  --protein \
  --total_length \
  --canonical  

conda activate vcf_bcf

#2) Extract relevant information
awk -F'\t' 'BEGIN {OFS="\t"}
NR > 1 {
    print $2, $3, $4, $5, $6, $7, $10, $11, $12, $16;
}' "${PATH_MSA_out}MSA_monomorphic_modified_vep.txt" > "${PATH_MSA_out}MSA_monomorphic_vep_modified_vep_final.txt"

#############################################################