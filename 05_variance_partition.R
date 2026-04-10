#############################################################
library("dplyr")
library("tidyverse")
library("MM4LMM")
#############################################################
# A) Paths, input and parameters
#############################################################

set.seed(0)

cores <- 25

submission <- "submission1"

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/variance_partition/")

data_list <- readRDS(paste0(path_input, "prepared_data/Phenotypes.rds"))
G_GRM <- readRDS(paste0(path_input, "kinships/random_effects/GRMs.rds"))[["Geno"]]
P_GRM_list <- readRDS(paste0(path_input, "kinships/random_effects/GRMs_DA.rds"))
P_GRM_permuted_list <- readRDS(paste0(path_input, "kinships/random_effects/GRMs_DA_permuted.rds"))
loads <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA.rds")) %>% distinct(id, G)

#############################################################
#B) Variance partition (unpermuted)
#############################################################

llik_output <- data.frame()
sigma_output <- data.frame()

for (current_trait in names(data_list)) {
  print(current_trait)

  data_subset <- data_list[[current_trait]] %>%
    dplyr::mutate_at(vars("id"), as.character) %>%
    drop_na(phenotype) %>%
    left_join(loads, by = "id")

  X <- as.formula("phenotype ~ 1 + G + PC1 + PC2 + PC3")

  Z_error <- diag(1, nrow(data_subset))

  G_GRM_ext <- G_GRM[data_subset$id, data_subset$id]

  M0_list <- list("G" = G_GRM_ext, "error" = Z_error)

  M0_meta <- list("trait" = current_trait, "permuted" = "no", "permutation" = "0", "q_interval" = "baseline", "model" = "M0")

  out_M0 <- MMEst(Y = data_subset$phenotype, Cofactor = data_subset, Method = "Reml", formula = X, VarList = M0_list, NbCores = cores)$NullModel

  M0_llik_output <- out_M0$`LogLik (Reml)` %>%
    as.data.frame() %>%
    rename(llik = ".") %>%
    dplyr::mutate(!!!M0_meta)

  llik_output <- rbind(M0_llik_output, llik_output)

  M0_sigma_output <- out_M0$Sigma2 %>%
    as.data.frame() %>%
    rename(sigma = ".") %>%
    rownames_to_column("component") %>%
    dplyr::mutate(!!!M0_meta)

  sigma_output <- rbind(M0_sigma_output, sigma_output)

  for (current_interval in names(P_GRM_list)) {
    P_GRM_ext <- P_GRM_list[[current_interval]][["P"]][data_subset$id, data_subset$id]

    M1_list <- list("G" = G_GRM_ext, "P" = P_GRM_ext, "error" = Z_error)

    M1_meta <- list("trait" = current_trait, "permuted" = "no", "permutation" = "0", "q_interval" = current_interval, "model" = "M1")

    tryCatch(
      {
        out_M1 <- MMEst(Y = data_subset$phenotype, Cofactor = data_subset, Method = "Reml", formula = X, VarList = M1_list, NbCores = cores)$NullModel

        llik_output <- bind_rows(llik_output, bind_rows(tibble(llik = out_M1$`LogLik (Reml)`)) %>% dplyr::mutate(!!!M1_meta))

        sigma_temp <- out_M1$Sigma2 %>%
          as.data.frame() %>%
          rownames_to_column("component") %>%
          setNames(c("component", "sigma")) %>%
          dplyr::mutate(!!!M1_meta)

        sigma_output <- rbind(sigma_output, sigma_temp)
      },
      error = function(e) {
        cat("Error:", e$message, "\n")
      }
    )
  }
}

saveRDS(llik_output, paste0(path_output, "Likelihoods_unpermuted_", submission, ".rds"))
saveRDS(sigma_output, paste0(path_output, "Variance_components_unpermuted_", submission, ".rds"))

#############################################################
#C) Variance partition (permuted)
#############################################################

llik_output <- data.frame()
sigma_output <- data.frame()

for (current_trait in names(data_list)) {
  print(current_trait)

  data_subset <- data_list[[current_trait]] %>%
    dplyr::mutate_at(vars("id"), as.character) %>%
    drop_na(phenotype) %>%
    left_join(loads, by = "id")

  X <- as.formula("phenotype ~ 1 + G + PC1 + PC2 + PC3")

  Z_error <- diag(1, nrow(data_subset))

  G_GRM_ext <- Matrix::nearPD(G_GRM[data_subset$id, data_subset$id])$mat %>% as.matrix()

  for (current_permute in names(P_GRM_permuted_list)) {
    P_GRM_permuted_list_upd <- P_GRM_permuted_list[[current_permute]]

    for (current_interval in names(P_GRM_permuted_list_upd)) {
      meta <- list("trait" = current_trait, permuted = "yes", "permutation" = current_permute, "q_interval" = current_interval, "model" = "M1")

      P_GRM_ext <- Matrix::nearPD(P_GRM_permuted_list_upd[[current_interval]][["P"]][data_subset$id, data_subset$id])$mat %>% as.matrix()

      M1_list <- list("G" = G_GRM_ext, "P" = P_GRM_ext, "error" = Z_error)

      tryCatch(
        {
          out_M1 <- MMEst(Y = data_subset$phenotype, Cofactor = data_subset, Method = "Reml", formula = X, VarList = M1_list, NbCores = cores)$NullModel

          llik_output <- bind_rows(llik_output, bind_rows(tibble(llik = out_M1$`LogLik (Reml)`)) %>% dplyr::mutate(!!!meta))

          sigma_output <- bind_rows(sigma_output, bind_rows(tibble(sigma = out_M1$Sigma2)) %>% rownames_to_column("component") %>% dplyr::mutate(!!!meta))
        },
        error = function(e) {
          cat("Error:", e$message, "\n")
        }
      )
    }
  }
}

saveRDS(llik_output, paste0(path_output, "Likelihoods_permuted_", submission, ".rds"))
saveRDS(sigma_output, paste0(path_output, "Variance_components_permuted_", submission, ".rds"))

#############################################################
