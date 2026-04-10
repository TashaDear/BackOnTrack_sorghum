#####################################
library("fastdfe")
fd <- load_fastdfe()
#####################################
# Paths
path_vcf <- "/home/tasha/backontrack/input/vcf_annotations/MSA/"
path_ref <- "/home/tasha/sorghum_vep/data/references/"
#####################################

# Annotate degeneracy
ann <- fd$Annotator(
  vcf = paste0(path_vcf, "MSA_merged_SAP_biallelic.vcf.gz"),
  fasta = paste0(path_ref, "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz"),
  gff = paste0(path_ref, "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"),
  annotations = c(fd$DegeneracyAnnotation()),
  output = paste0(path_vcf, "MSA_merged_SAP_biallelic_deg.vcf.gz")
)

fd$Annotator$annotate(ann)

# Ancestral allele annotation -> outgroups specified based on increasing divergence
ann <- fd$Annotator(
  vcf = paste0(path_vcf, "MSA_merged_SAP_biallelic_deg.vcf.gz"),
  annotations = c(fd$MaximumLikelihoodAncestralAnnotation(outgroups = c("Erianthus_rufipilus", "Zea_mays_1"), n_ingroups = 20)),
  output = paste0(path_vcf, "MSA_merged_SAP_biallelic_deg_AA.vcf.gz")
)

fd$Annotator$annotate(ann)
#####################################