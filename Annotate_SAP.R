#####################################
library("fastdfe")
fd <- load_fastdfe()
#####################################
# Paths

path_vcf <- "/home/tasha/backontrack/input/vcf_annotations/SAP/"
path_ref <- "/home/tasha/sorghum_vep/data/references/"

#####################################

# Annotate degeneracy
ann <- fd$Annotator(
  vcf = paste0(path_vcf, "SAP.imputed_upd_CDS.vcf.gz"),
  fasta = paste0(path_ref, "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz"),
  gff = paste0(path_ref, "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"),
  annotations = c(fd$DegeneracyAnnotation()),
  output = paste0(path_vcf, "SAP.imputed_upd_CDS_deg.vcf.gz")
)

fd$Annotator$annotate(ann)
#####################################