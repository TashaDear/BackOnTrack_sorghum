#############################################################
library("dplyr")
library("tidyverse")
library("MM4LMM")
#############################################################
# A) Paths
#############################################################

set.seed(0)

cores <- 25

submission <- "submission1"

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/mean_partition/")

data_list <- readRDS(paste0(path_input, "prepared_data/Phenotypes.rds"))
G_GRM <- readRDS(paste0(path_input, "kinships/random_effects/GRMs.rds"))[["Geno"]]
loads <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA.rds"))
loads_perm <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA_permuted.rds"))
#############################################################
# B) Mean partition (unpermuted)
#############################################################

final_output <- data.frame()

for (current_trait in names(data_list)) {
  print(current_trait)

  data_subset <- data_list[[current_trait]] %>%
    mutate(id = as.character(id)) %>%
    drop_na(phenotype)

  meta <- list("trait" = current_trait, "permuted" = "no", "permutation" = "0")

  for (current_interval in unique(loads$esm_interval)) {
    print(current_interval)

    loads_final <- loads %>% filter(esm_interval == current_interval)

    meta_upd <- c(meta, list("q_interval" = unique(loads_final$q_interval), "esm_interval" = current_interval))

    data_updated <- data_subset %>% left_join(loads_final, by = "id")

    V_list <- list("G" = G_GRM[data_updated$id, data_updated$id], "error" = diag(1, nrow(data_updated)))

    out <- MMEst(
      Y = data_updated$phenotype, Cofactor = data_updated,
      formula = as.formula("phenotype ~ 1 + G + P + PC1 + PC2 + PC3"), Method = "ML", VarList = V_list, NbCores = cores
    )

    Beta <- out$NullModel$Beta
    Sigma <- out$NullModel$VarBeta
    colnames(Sigma) <- rownames(Sigma)
    sig_P <- aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) == "P"))$result$chi2["P"]

    temp_output <- data.frame("sig" = sig_P, "effect" = Beta["P"], "se" = sqrt(Sigma["P", "P"])) %>% mutate(!!!meta_upd)

    final_output <- bind_rows(final_output, temp_output)
  }
}

saveRDS(final_output, paste0(path_output, "Mean_partition_", submission, ".rds"))

#############################################################
# C) Mean partition (permuted)
#############################################################

final_output <- data.frame()

for (current_trait in names(data_list)) {
  print(current_trait)

  data_subset <- data_list[[current_trait]] %>%
    mutate(id = as.character(id)) %>%
    drop_na(phenotype)

  X_formula <- as.formula("phenotype ~ 1 + G + P + PC1 + PC2 + PC3")

  for (current_interval in unique(loads_perm$esm_interval)) {
    V_list <- list("G" = G_GRM[data_subset$id, data_subset$id], "error" = diag(1, nrow(data_subset)))

    loads_perm_upd <- loads_perm %>% dplyr::filter(esm_interval == current_interval)

    perm_list <- loads_perm_upd %>%
      dplyr::select(c("id", "G", "P", "permutation_seed")) %>%
      split(.$permutation_seed) %>%
      lapply(function(input) {
        rownames(input) <- NULL
        mat <- input %>%
          column_to_rownames("id") %>%
          dplyr::select(-permutation_seed) %>%
          as.matrix()
        mat[data_subset$id, , drop = FALSE]
      })

    result <- tryCatch(
      {
        out <- MMEst(
          Y = data_subset$phenotype, Cofactor = data_subset, X = perm_list,
          formula = X_formula, Method = "ML", VarList = V_list, NbCores = cores
        )

        lapply(names(out), function(current_permute) {
          Beta <- out[[current_permute]]$Beta
          Sigma <- out[[current_permute]]$VarBeta
          colnames(Sigma) <- rownames(Sigma)

          sig_P <- aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) %in% c("P")))$result$chi2["P"]

          data.frame(
            "trait" = current_trait, "sig" = sig_P, "effect" = Beta["P"], "se" = sqrt(Sigma["P", "P"]),
            "permuted" = "yes", "permutation" = current_permute,
            "q_interval" = unique(loads_perm_upd$q_interval), "esm_interval" = current_interval
          )
        })
      },
      error = function(e) {
        cat("Error:", e$message, "\n")
        return(NULL)
      }
    )

    final_output <- bind_rows(final_output, result)
  }
}

saveRDS(final_output, paste0(path_output, "Mean_partition_permuted_", submission, ".rds"))

#############################################################
