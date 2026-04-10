#############################################################
library(AGHmatrix)
library(ASRgenomics)
library(caret)
library(data.table)
library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(readr)
library(rlang)
library(SNPRelate)
library(tidyverse)
library(viridis)
library(purrr)
#############################################################
# A) Paths
#############################################################

set.seed(0)

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/")
path_results <- paste0(path, "results/")

vcf_in <- paste0(path_input, "SAP/SAP_formatted/SAP.imputed_upd.vcf.gz")
GDS_out <- paste0(path_input, "SAP/SAP_formatted/SAP.imputed_upd_biallelic.gds")

permutation_seeds <- 1:5000

#############################################################
# B) Load input
#############################################################

passport <- read_csv(paste0(path_input, "SAP/tpj15853-sup-0001-files1.csv")) %>%
  as.data.frame() %>%
  rename("id" = "Taxa") %>%
  filter(!is.na(K.Cluster))

Missense_list <- readRDS(paste0(path_input, "/Missense_list.rds"))

SNPs_AA_missense <- Missense_list[["Missense_AA_df"]]
SNPs_missense <- Missense_list[["Missense_AA_df"]] %>% dplyr::filter(Degeneracy %in% c("0"), !is.na(esm_score))

#############################################################
# C) Create functions
#############################################################
# Function #1
get_geno <- function(geno_file, snp.ids, SNPs_position) {
  geno <- snpgdsGetGeno(geno_file, snp.id = snp.ids, .snpread = NA, with.id = TRUE, verbose = FALSE)
  rownames(geno$genotype) <- geno$sample.id
  colnames(geno$genotype) <- geno$snp.id
  geno_ALT <- 2 - geno$genotype
  colnames(geno_ALT) <- SNPs_position$CHROM_POS[match(colnames(geno_ALT), as.character(SNPs_position$snp.id))]
  return(geno_ALT)
}

#############################################################
# Function #2
create_GRM <- function(geno_list) {
  GRM_list <- list()
  scaling_df <- data.frame()

  for (geno_name in names(geno_list)) {
    sum_2pq <- sum(2 * (colMeans(geno_list[[geno_name]]) / 2) * (1 - (colMeans(geno_list[[geno_name]]) / 2)))
    X_center <- geno_list[[geno_name]] %>% scale(., center = T, scale = F)
    GRM_list[[geno_name]] <- tcrossprod(X_center) / sum_2pq
    scaling_df <- rbind(scaling_df, data.frame("GRM" = geno_name, "Scaling" = sum_2pq))
  }
  GRM_list[["Scaling"]] <- scaling_df
  return(GRM_list)
}

#############################################################
# Function 3
create_GRM_prioritized <- function(input, geno_matrix, polarization, permute, seed, quantiles) {
  set.seed(seed)

  GRM_list <- list()

  pol_map <- c("BA" = "esm_score_BA", "DA" = "esm_score_DA")

  input_upd <- input %>%
    dplyr::mutate(
      POL = case_when(
        polarization == "BA" ~ if_else(esm_score >= 0, ALT, "REF"),
        polarization == "DA" ~ if_else(ALT != AA, ALT, "REF"), TRUE ~ NA_character_
      ),
      esm_score_pol = .data[[pol_map[polarization]]]
    )

  input_ALT_is_DA <- input_upd %>%
    dplyr::filter(ALT != AA) %>%
    dplyr::mutate(esm_score_perm = sample(esm_score, replace = FALSE))

  input_REF_is_DA <- input_upd %>%
    dplyr::filter(ALT == AA) %>%
    dplyr::mutate(esm_score_perm = sample(esm_score, replace = FALSE))

  input_final <- rbind(input_ALT_is_DA, input_REF_is_DA) %>%
    dplyr::mutate(
      esm_score_perm_pol = esm_score_perm * ifelse(POL == ALT, 1, -1),
      score = if (permute == "yes") esm_score_perm_pol else esm_score_pol,
      score_norm = score / max(abs(score))
    )

  for (i in 1:(length(quantiles) - 1)) {
    q_interval <- paste0("[", names(quantiles)[i], ", ", names(quantiles)[i + 1], ")")

    sites <- input_final %>%
      dplyr::filter(score >= quantiles[i], if (i < length(quantiles) - 1) score < quantiles[i + 1] else score <= quantiles[i + 1])

    P_list <- list("P" = geno_matrix[, sites$CHROM_POS])

    GRM_list[[q_interval]] <- create_GRM(P_list)
  }

  return(GRM_list)
}

### Helper function
process_GRMs <- function(input, geno_matrix, polarization, quantiles) {
  GRM_list <- create_GRM_prioritized(input, geno_matrix, polarization, "no", 0, quantiles)
  saveRDS(GRM_list, paste0(path_input, "kinships/random_effects/GRMs_", polarization, ".rds"))

  GRM_list_perm <- map(permutation_seeds, ~ create_GRM_prioritized(input, geno_matrix, polarization, "yes", .x, quantiles)) %>%
    set_names(paste0("perm_", permutation_seeds))
  saveRDS(GRM_list_perm, paste0(path_input, "kinships/random_effects/GRMs_", polarization, "_permuted.rds"))
}

#############################################################
# Function 4
compute_loads <- function(input, geno_matrix, polarization, permute, seed, quantiles) {
  set.seed(seed)

  output <- data.frame()

  pol_map <- c("BA" = "esm_score_BA", "DA" = "esm_score_DA")

  input_upd <- input %>%
    dplyr::mutate(
      POL = case_when(
        polarization == "BA" ~ if_else(esm_score >= 0, ALT, "REF"),
        polarization == "DA" ~ if_else(ALT != AA, ALT, "REF"), TRUE ~ NA_character_
      ),
      esm_score_pol = .data[[pol_map[polarization]]]
    )

  input_ALT_is_DA <- input_upd %>%
    dplyr::filter(ALT != AA) %>%
    dplyr::mutate(esm_score_perm = sample(esm_score, replace = FALSE))

  input_REF_is_DA <- input_upd %>%
    dplyr::filter(ALT == AA) %>%
    dplyr::mutate(esm_score_perm = sample(esm_score, replace = FALSE))

  input_final <- rbind(input_ALT_is_DA, input_REF_is_DA) %>%
    dplyr::mutate(
      esm_score_perm_pol = esm_score_perm * ifelse(POL == ALT, 1, -1),
      score = if (permute == "yes") esm_score_perm_pol else esm_score_pol,
      score_norm = score / max(abs(score))
    )

  geno_full <- geno_matrix[, input_final$CHROM_POS] %>% as.data.frame()

  recode <- input_final %>%
    dplyr::filter(POL != ALT) %>%
    pull(CHROM_POS)

  geno_recoded <- geno_full

  geno_recoded[, recode] <- 2 - geno_full[, recode]

  load_G <- geno_recoded[, ] %>%
    {
      rowSums(.) %>% as.data.frame()
    } %>%
    rownames_to_column("id") %>%
    rename("G" = ".")

  if (permute == "no") {
    weights_vec <- input_final$score_norm[match(colnames(geno_recoded), input_final$CHROM_POS)]

    load_total <- sweep(geno_recoded, 2, weights_vec, `*`) %>%
      {
        rowSums(.) %>% as.data.frame()
      } %>%
      rownames_to_column("id") %>%
      rename("P" = ".")

    saveRDS(load_total, paste0(path_input, "kinships/fixed_effects/GW_load_DA_weighted.rds"))
  }

  for (i in 1:(length(quantiles) - 1)) {
    esm_interval <- paste0("[", round(quantiles[i], 1), ", ", round(quantiles[i + 1], 1), ")")
    q_interval <- paste0("[", as.numeric(names(quantiles))[i], ", ", as.numeric(names(quantiles))[i + 1], ")")

    sites <- input_final %>%
      dplyr::filter(score >= quantiles[i], if (i < length(quantiles) - 1) score < quantiles[i + 1] else score <= quantiles[i + 1])

    load_P <- geno_recoded[, sites$CHROM_POS] %>%
      {
        rowSums(.) %>% as.data.frame()
      } %>%
      rownames_to_column("id") %>%
      rename("P" = ".") %>%
      dplyr::mutate(
        "esm_interval" = esm_interval, "q_interval" = q_interval, "sites" = length(unique(sites$CHROM_POS)),
        "polarization" = polarization, "permute" = permute, "permutation_seed" = paste0("perm_", seed)
      )

    output <- bind_rows(load_G %>% left_join(load_P, by = "id"), output)
  }

  return(output)
}

# Helper function
process_loads <- function(input, geno_matrix, polarization, quantiles) {
  GW_loads <- compute_loads(input, geno_matrix, polarization, "no", 0, quantiles)
  saveRDS(GW_loads, paste0(path_input, "kinships/fixed_effects/GW_load_", polarization, ".rds"))

  GW_loads_perm <- purrr::map_dfr(permutation_seeds, ~ compute_loads(input, geno_matrix, polarization, "yes", .x, quantiles))
  saveRDS(GW_loads_perm, paste0(path_input, "kinships/fixed_effects/GW_load_", polarization, "_permuted.rds"))
}

#############################################################
# Function 5
get_named_quantiles <- function(data, column, probs) {
  q <- quantile(data[[column]], probs = probs, na.rm = TRUE)

  setNames(as.numeric(q), probs)
}

#############################################################
# D) Load VCF
#############################################################

snpgdsVCF2GDS(vcf_in, GDS_out, method = "biallelic.only")

geno_file <- snpgdsOpen(GDS_out)

SNPs_position <- snpgdsSNPList(geno_file) %>%
  as.data.frame() %>%
  dplyr::mutate(CHROM_POS = paste(chromosome, position, sep = "_"))
SNPs_position_MIS <- SNPs_position %>% dplyr::filter(CHROM_POS %in% SNPs_missense$CHROM_POS)

Geno <- get_geno(geno_file, snp.id = NULL, SNPs_position = SNPs_position)
Geno_MIS <- get_geno(geno_file, snp.id = SNPs_position_MIS$snp.id, SNPs_position = SNPs_position)

prioritized_SNPs <- SNPs_position_MIS %>%
  dplyr::filter(CHROM_POS %in% colnames(Geno_MIS)) %>%
  left_join(SNPs_missense, by = "CHROM_POS") %>%
  dplyr::mutate(
    snp.id = as.character(snp.id),
    afreq_REF = afreq,
    afreq_ALT = (1 - afreq_REF),
    esm_score_BA = abs(esm_score)
  )

#############################################################
# E) Prepare data frames
############################################################

AA_df <- prioritized_SNPs %>%
  dplyr::filter(Degeneracy %in% c("0"), !is.na(AA), !is.na(esm_score)) %>%
  dplyr::mutate(esm_score_DA = if_else(ALT == AA, -esm_score, esm_score))

SNPs_AA_missense_upd <- SNPs_AA_missense %>% left_join(SNPs_position, by = "CHROM_POS")

#############################################################
# F) Prepare quantiles
#############################################################

current_probs <- seq(0, 1, by = 0.1)

quantiles_list <- tibble(quantile = current_probs) %>%
  bind_cols(
    "DA" = get_named_quantiles(AA_df, "esm_score_DA", current_probs),
    "BA" = get_named_quantiles(AA_df, "esm_score_BA", current_probs)
  )
saveRDS(quantiles_list, paste0(path_input, "kinships/fixed_effects/Quantiles_df_final.rds"))

AA_df_upd <- AA_df %>%
  dplyr::mutate(esm_interval_DA = cut(esm_score_DA,
    breaks = quantiles_list[["DA"]],
    labels = paste0(
      "[", head(round(quantiles_list[["DA"]], 1), -1), ", ",
      tail(round(quantiles_list[["DA"]], 1), -1), ")"
    ), include.lowest = TRUE, right = FALSE
  )) %>%
  dplyr::mutate(esm_interval_BA = cut(esm_score_BA,
    breaks = quantiles_list[["BA"]],
    labels = paste0(
      "[", head(round(quantiles_list[["BA"]], 1), -1), ", ",
      tail(round(quantiles_list[["BA"]], 1), -1), ")"
    ), include.lowest = TRUE, right = FALSE
  )) %>%
  dplyr::select(-c("position", "allele", "GENE_ID", "TRANSCRIPT_ID", "afreq_REF", "REF", "snp.id", "chromosome"))
saveRDS(AA_df_upd, paste0(path_input, "kinships/fixed_effects/Score_intervals_final.rds"))

#############################################################
# G) Prepare random and fixed effects
#############################################################

geno_list <- list("Geno" = Geno, "Missense" = Geno_MIS)
saveRDS(geno_list, paste0(path_input, "kinships/Geno_list.rds"))

GRM_list <- create_GRM(geno_list)
saveRDS(GRM_list, paste0(path_input, "kinships/random_effects/GRMs.rds"))

# Random effects
process_GRMs(AA_df, Geno_MIS, "DA", quantiles_list[["DA"]])

# Fixed effects
process_loads(AA_df, Geno_MIS, "DA", quantiles_list[["DA"]])

#############################################################
# I) Determine SNPs with MAF above lower threshold
#############################################################

set.seed(0)

MAF_across <- 0.05
MAF_within <- 0.01

SNPs_position_across <- SNPs_position %>%
  dplyr::mutate(REF_freq = afreq, ALT_freq = 1 - afreq, chromosome = as.character(chromosome)) %>%
  dplyr::filter(ALT_freq > MAF_across, REF_freq > MAF_across, CHROM_POS %in% colnames(Geno))

sample_clusters <- setNames(passport$K.Cluster, passport$id)
valid_samples <- intersect(rownames(Geno), names(sample_clusters))
Geno_sub <- Geno[valid_samples, SNPs_position_across$snp.id, drop = FALSE]

sample_clusters_dt <- data.table(ID = names(sample_clusters), K.Cluster = sample_clusters)
geno_dt <- data.table(ID = rownames(Geno_sub), Geno_sub)
geno_dt <- merge(geno_dt, sample_clusters_dt, by = "ID")
setkey(geno_dt, K.Cluster)

AF_dt <- geno_dt[, lapply(.SD, \(x) mean(x, na.rm = TRUE) / 2), by = K.Cluster, .SDcols = setdiff(names(geno_dt), c("ID", "K.Cluster"))]
MAF_pr_cluster <- melt(AF_dt, id.vars = "K.Cluster", variable.name = "CHROM_POS", value.name = "AF")

SNPs_to_keep <- MAF_pr_cluster %>%
  group_by(CHROM_POS) %>%
  dplyr::filter(all(AF > MAF_within)) %>%
  ungroup() %>%
  pull(CHROM_POS)

SNPs_position_within <- SNPs_position_across %>% dplyr::filter(CHROM_POS %in% SNPs_to_keep)

input_LD_list <- list(
  "Geno_sub" = Geno[, unique(c(AA_df$CHROM_POS, SNPs_position_across$CHROM_POS))],
  "SNPs_across" = SNPs_position_across,
  "SNPs_within" = SNPs_position_within
)

saveRDS(input_LD_list, paste0(path_input, "LDdecay/Input_LD_list.rds"))

#############################################################
# J) Determine PCs and prepare dataframe
#############################################################

PCs_geno <- prcomp(Geno, center = T, scale = F)
PCs_geno_upd <- PCs_geno[c("x", "sdev")]
saveRDS(PCs_geno_upd, paste0(path_input, "prepared_data/PCA_data.rds"))

PCs_unscaled <- PCs_geno$x %>%
  as.data.frame() %>%
  dplyr::select(c("PC1", "PC2", "PC3")) %>%
  rownames_to_column("id")

passport <- read_csv(paste0(path_input, "SAP/tpj15853-sup-0001-files1.csv")) %>%
  as.data.frame() %>%
  dplyr::rename("id" = "Taxa")

phenotypes <- read_csv(paste0(path_input, "SAP/SAP.pheno.csv")) %>%
  as.data.frame() %>%
  dplyr::rename("id" = "ID") %>%
  left_join(passport, by = "id") %>%
  left_join(PCs_unscaled, by = "id") %>%
  pivot_longer(
    cols = c("Amylose", "Fat", "Cal.g", "Starch", "Protein", "DTA", "PH", "GN", "GW", "GY", "FLH", "PL", "BL"),
    names_to = "trait", values_to = "phenotype"
  ) %>%
  dplyr::filter(id %in% rownames(Geno) & !is.na(K.Cluster)) %>%
  dplyr::mutate(
    K.Cluster = as.character(K.Cluster),
    id = as.character(id)
  ) %>%
  split(.$trait)
saveRDS(phenotypes, paste0(path_input, "prepared_data/Phenotypes.rds"))

#############################################################
snpgdsClose(geno_file)
#############################################################
