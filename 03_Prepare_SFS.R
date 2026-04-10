#############################################################
library("dplyr")
library("ggplot2")
library("readr")
library("SNPRelate")
library("tidyverse")
library("viridis")
library("fastdfe")
library("data.table")
library("foreach")
library("reticulate")
library("ggridges")
library("doParallel")
registerDoParallel(cores = detectCores() - 1)
fd <- load_fastdfe()
np <- import("numpy")
#############################################################
# A) Paths
#############################################################

set.seed(0)

submission <- "submission1"

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/")
path_results <- paste0(path, "results/")

vcf_in <- paste0(path_input, "SAP/SAP_formatted/SAP.imputed_upd.vcf.gz")
GDS_out <- paste0(path_input, "SAP/SAP_formatted/SAP.imputed_upd_biallelic.gds")

#############################################################
# B) Format the input for analysis
#############################################################

# Missense sites
Missense_list <- readRDS(paste0(path_input, "/Missense_list.rds"))

# Evolutionary scores
esm_scores_potential_df <- fread(paste0(path_input, "esm_predictions/NJ_pred/Sorghum_esm_scores.csv"),
  col.names = c("TRANSCRIPT_ID", "mutation", "esm_score")
)[
  ,
  `:=`(
    REF_aminoacid = substr(mutation, 1, 1),
    aminoacid_position = substr(mutation, 2, nchar(mutation) - 1),
    ALT_aminoacid = substr(mutation, nchar(mutation), nchar(mutation))
  )
][, .(TRANSCRIPT_ID, mutation, REF_aminoacid, aminoacid_position, ALT_aminoacid, esm_score)]

# Monomorphic sites (alignable to outgroups)
Mono_sites <- fread(paste0(path_input, "MSA/MSA_formatted/MSA_assumed_monomorphic_modified_vep.txt"), sep = "\t", header = TRUE) %>%
  setDT() %>%
  .[, c("CHROM", "POS") := tstrsplit(Location, ":")] %>%
  .[, `:=`(
    CHROM_POS = paste(CHROM, POS, sep = "_"),
    REF_codon = sub("/.*", "", Codons),
    REF_aminoacid = substr(Amino_acids, 1, 1),
    Protein_position = sub("/.*", "", Protein_position),
    effect = sub("_.*", "", Consequence)
  )] %>%
  dplyr::filter(effect != "synonymous") %>%
  .[nchar(Codons) >= 3] %>%
  dplyr::select(c("CHROM_POS", "CHROM", "POS", "Feature", "REF_codon", "REF_aminoacid", "Protein_position")) %>%
  rename("TRANSCRIPT_ID" = "Feature", "aminoacid_position" = "Protein_position") %>%
  left_join(esm_scores_potential_df, by = c("TRANSCRIPT_ID", "aminoacid_position", "REF_aminoacid")) %>%
  filter(!is.na(esm_score)) %>%
  setDT() %>%
  {
    set(., j = "Identifier", value = paste(.$TRANSCRIPT_ID, .$CHROM_POS, sep = "_"))
    .
  } %>%
  dplyr::select(-c("mutation", "aminoacid_position"))

Mono_sites <- Mono_sites %>%
  group_by(CHROM_POS) %>%
  dplyr::filter(n_distinct(TRANSCRIPT_ID) == 1) %>%
  ungroup() %>%
  split(.$CHROM)

#
Scores <- readRDS(paste0(path_input, "kinships/fixed_effects/Score_intervals_final.rds")) %>%
  dplyr::select(CHROM_POS, AA, ALT, esm_score, esm_score_DA, esm_interval_DA) %>%
  dplyr::mutate(esm_interval_DA = recode(esm_interval_DA, "[3.9, 11.3)" = "[3.9, 11.3]")) %>%
  dplyr::filter(esm_score_DA != "[-12.8, 11.3]")

#############################################################
# C) Determine target sites
#############################################################

out <- fd$Parser(
  n = 10,
  max_sites = 100000,
  vcf = paste0(path_input, "MSA/MSA_formatted/Ingroups_deg_AA_trimmed_biallelic.vcf.gz"),
  fasta = paste0(
    path_input,
    "references/version_3/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz"
  ),
  gff = paste0(
    path_input,
    "references/version_3/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3.gz"
  ),
  filtrations = c(fd$SNPFiltration()),
  annotations = c(fd$DegeneracyAnnotation()),
  stratifications = c(fd$DegeneracyStratification()),
  target_site_counter = fd$TargetSiteCounter(n_samples = 10000, n_target_sites = 1000000)
)

sfs <- fd$Parser$parse(out)
Ratio <- sfs$data["neutral"][1, ] / sfs$data["selected"][1, ]

# Polymorphic in population, 0fold
P_NS <- Missense_list[["Missense_df"]] %>%
  filter(Degeneracy %in% c("0")) %>%
  distinct(CHROM_POS) %>%
  count() %>%
  pull()

Target_site_count_sel <- P_NS * 1000 # increased from 1000 to 10000 04/02/26 in the upd file

#############################################################
# D) functions
#############################################################

# Codons and aminoacids table
Codons_DNA <- tibble::tribble(
  ~codon, ~amino_acid,
  "TTT", "F", "TTC", "F", "TTA", "L", "TTG", "L",
  "CTT", "L", "CTC", "L", "CTA", "L", "CTG", "L",
  "ATT", "I", "ATC", "I", "ATA", "I", "ATG", "M",
  "GTT", "V", "GTC", "V", "GTA", "V", "GTG", "V",
  "TCT", "S", "TCC", "S", "TCA", "S", "TCG", "S",
  "CCT", "P", "CCC", "P", "CCA", "P", "CCG", "P",
  "ACT", "T", "ACC", "T", "ACA", "T", "ACG", "T",
  "GCT", "A", "GCC", "A", "GCA", "A", "GCG", "A",
  "TAT", "Y", "TAC", "Y", "TAA", "Stop", "TAG", "Stop",
  "CAT", "H", "CAC", "H", "CAA", "Q", "CAG", "Q",
  "AAT", "N", "AAC", "N", "AAA", "K", "AAG", "K",
  "GAT", "D", "GAC", "D", "GAA", "E", "GAG", "E",
  "TGT", "C", "TGC", "C", "TGA", "Stop", "TGG", "W",
  "CGT", "R", "CGC", "R", "CGA", "R", "CGG", "R",
  "AGT", "S", "AGC", "S", "AGA", "R", "AGG", "R",
  "GGT", "G", "GGC", "G", "GGA", "G", "GGG", "G"
)

# Generate possible aminoacids based on REF codon
generate_codons <- function(codon) {
  codon_split <- strsplit(codon, NULL)[[1]]

  position <- which(codon_split == toupper(codon_split))

  possible_codons <- paste0(codon_split[1:(position - 1)], c("A", "C", "G", "T"), codon_split[(position + 1):length(codon_split)])

  amino_acids <- Codons_DNA %>%
    filter(codon %in% unique(toupper(possible_codons))) %>%
    pull(amino_acid)

  return(amino_acids)
}

# Filter out aminoacids that cannot occur with single-step mutations
Site_mutations <- function(input) {
  output <- setDT(input)[,
    {
      ALT_aminoacid <- unlist(generate_codons(REF_codon))
      .(ALT_aminoacid)
    },
    by = Identifier
  ][, .SD[n_distinct(ALT_aminoacid) >= 3], by = Identifier][ALT_aminoacid != "Stop"]

  final_output <- input[output, on = .(Identifier, ALT_aminoacid), nomatch = 0][REF_aminoacid != ALT_aminoacid]

  return(final_output)
}

# Determine esm scores for monomorphic sites
Monomorphic_esm_mean <- function(input) {
  Monomorphic_list <- list()

  for (Chr in names(input)) {
    print(paste("Chromosome:", Chr))

    Monomorphic_list[[Chr]] <- Site_mutations(input[[Chr]])[,
      .(esm_score_mean = sum(esm_score, na.rm = TRUE) / .N),
      by = .(CHROM_POS, TRANSCRIPT_ID)
    ][!is.na(esm_score_mean)]
  }

  Monomorphic_df <- Monomorphic_list %>% do.call(rbind, .)

  saveRDS(Monomorphic_df, paste0(path_input, "prepared_data/Monomorphic_sites_esm_scores_mean_df.rds"))
}

#############################################################
# E) Run functions
#############################################################

Monomorphic_esm_mean(Mono_sites)

#############################################################
# F) Prepare geno input
#############################################################

snpgdsVCF2GDS(vcf_in, GDS_out, method = "biallelic.only")

geno_file <- snpgdsOpen(GDS_out)

SNPs_position <- snpgdsSNPList(geno_file) %>%
  as.data.frame() %>%
  dplyr::mutate(CHROM_POS = paste(chromosome, position, sep = "_"))

snps <- SNPs_position %>%
  dplyr::filter(CHROM_POS %in% Missense_list[["Missense_AA_df"]]$CHROM_POS) %>%
  pull(snp.id)

Geno <- snpgdsGetGeno(geno_file, snp.id = snps, .snpread = NA, with.id = T, verbose = F)

rownames(Geno$genotype) <- Geno$sample.id
colnames(Geno$genotype) <- Geno$snp.id
colnames(Geno$genotype) <- SNPs_position$CHROM_POS[match(colnames(Geno$genotype), SNPs_position$snp.id)]

SNPs_recode <- Missense_list[["Missense_AA_df"]] %>%
  dplyr::filter(ALT == AA) %>%
  pull(CHROM_POS) %>%
  as.character()

Geno_ALT <- 2 - Geno$genotype

Geno_DA <- Geno_ALT

Geno_DA[, SNPs_recode] <- 2 - Geno_ALT[, SNPs_recode]

#############################################################
# G) Prepare input for SFS
#############################################################

Allele_freq <- Geno_DA %>%
  colSums() %>%
  as.data.frame() %>%
  rename("Derived_count" = ".") %>%
  dplyr::mutate("AF" = Derived_count / (2 * nrow(Geno_DA))) %>%
  rownames_to_column("CHROM_POS") %>%
  left_join(Missense_list[["Missense_AA_df"]], by = c("CHROM_POS")) %>%
  dplyr::mutate(esm_score_polarized = if_else(ALT == AA, -esm_score, esm_score))

SFS_input <- Allele_freq %>%
  dplyr::filter(!is.na(AA)) %>%
  dplyr::mutate("AF_bin" = cut(AF,
    breaks = seq(0.00, 1, by = 0.05),
    labels = sprintf(
      "%.2f - %.2f", head(seq(0, 1, by = 0.05), -1),
      tail(seq(0, 1, by = 0.05), -1)
    ), include.lowest = TRUE
  ))

sites_4fold <- SFS_input %>% dplyr::filter(Degeneracy %in% c("4"))

sites_0fold <- SFS_input %>% dplyr::filter(!is.na(esm_score), Degeneracy %in% c("0"))

quantiles_vec <- quantile(sites_0fold$esm_score_polarized, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
setNames(as.numeric(quantiles_vec), seq(0, 1, by = 0.1))

SFS_input <- sites_0fold %>%
  dplyr::mutate(esm_range = cut(esm_score_polarized,
    breaks = quantiles_vec,
    labels = paste0(
      "[", head(round(quantiles_vec, 1), -1), ", ",
      tail(round(quantiles_vec, 1), -1), ")"
    ), include.lowest = TRUE, right = FALSE
  )) %>%
  split(.$esm_range)

fullrange <- paste0("[", round(quantiles_vec[1], 1), ", ", round(quantiles_vec[11], 1), "]")

SFS_input[[fullrange]] <- sites_0fold %>%
  dplyr::mutate(esm_range = paste0("[", round(quantiles_vec[1], 1), ", ", round(quantiles_vec[10], 1), "]"))

#############################################################
# H) Make SFS stratefied on evolutionary scores
#############################################################

prepare_SFS <- function(mean_or_random) {
  if (mean_or_random == "mean") {
    Mono_esm_distribution <- readRDS(paste0(path_input, "prepared_data/Monomorphic_sites_esm_scores_mean_df.rds")) %>%
      dplyr::rename(esm_score = esm_score_mean)
  } else if (mean_or_random == "random") {
    Mono_esm_distribution <- readRDS(paste0(path_input, "prepared_data/Monomorphic_sites_esm_scores_random_df.rds"))
  }

  total_monomorphic_sites <- length(unique(Mono_esm_distribution$CHROM_POS))

  SFS_list <- list()

  for (current_range in names(SFS_input)) {
    monomorphic_sites <- Mono_esm_distribution %>%
      dplyr::filter(esm_score >= min(SFS_input[[current_range]]$esm_score_polarized) &
        esm_score < max(SFS_input[[current_range]]$esm_score_polarized)) %>%
      pull(CHROM_POS) %>%
      length()

    Target_site_count_sel_upd <- (monomorphic_sites / total_monomorphic_sites) * Target_site_count_sel

    Target_site_count_neu_upd <- Ratio * Target_site_count_sel_upd

    SFS_sel <- SFS_input[[current_range]] %>%
      group_by(AF_bin) %>%
      summarise(Count = n(), .groups = "drop") %>%
      complete(AF_bin, fill = list(Count = 0)) %>%
      arrange(AF_bin) %>%
      pull(Count)

    SFS_neu <- sites_4fold %>%
      group_by(AF_bin) %>%
      summarise(Count = n(), .groups = "drop") %>%
      arrange(AF_bin) %>%
      pull(Count)

    SFS_list[[current_range]][["sel"]] <- c(round(Target_site_count_sel_upd), SFS_sel, 0)

    SFS_list[[current_range]][["neu"]] <- c(round(Target_site_count_neu_upd), SFS_neu, 0)
  }

  return(SFS_list)
}

#############################################################

SFS_mean_list <- prepare_SFS("mean")
saveRDS(SFS_mean_list, paste0(path_output, "SFS/SFS_list_input_mean_upd_", submission, ".rds"))

#############################################################
# I) Make csv file
#############################################################

SFS_list <- list("mean" = SFS_mean_list)

for (mean_or_random in names(SFS_list)) {
  input <- SFS_list[[mean_or_random]]

  df_list <- list()

  for (current_range in names(input)) {
    sel_vec <- input[[current_range]][["sel"]]

    neu_vec <- input[[current_range]][["neu"]]

    df_list[[current_range]] <- data.frame(esm_range = current_range, rbind("Selected" = sel_vec, "Neutral" = neu_vec), check.names = FALSE) %>%
      rownames_to_column("spectra")
  }

  final_df <- bind_rows(df_list)

  colnames(final_df) <- c(
    "SFS", "esm_range", "Monomorphic_ancestral",
    "[0.00-0.05)", "[0.05-0.10)", "[0.10-0.15)", "[0.15-0.20)", "[0.20-0.25)",
    "[0.25-0.30)", "[0.30-0.35)", "[0.35-0.40)", "[0.40-0.45)", "[0.45-0.50)",
    "[0.50-0.55)", "[0.55-0.60)", "[0.60-0.65)", "[0.65-0.70)", "[0.70-0.75)",
    "[0.75-0.80)", "[0.80-0.85)", "[0.85-0.90)", "[0.90-0.95)", "[0.95, 1]", "Monomorphic_derived"
  )

  write.csv(final_df, paste0(path_output, "SFS/SFS_", mean_or_random, "_input_upd.csv"), row.names = FALSE)
}

#############################################################
#############################################################
