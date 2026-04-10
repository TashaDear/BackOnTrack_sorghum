#############################################################
library("caret")
library("dplyr")
library("foreach")
library("MM4LMM")
library("tibble")
library("tidyr")
library("purrr")
library("qgg")
#############################################################
# A) Paths, input and parameters
#############################################################

set.seed(0)

cores <- 25

submission <- "submission1"

path_input <- "/home/tasha/backontrack/input/"
path_output <- "/home/tasha/backontrack/output/"

data_list <- readRDS(file.path(path_input, "prepared_data", "Phenotypes.rds")) %>%
  do.call(rbind, .) %>%
  dplyr::mutate(genetic_cluster = as.character(K.Cluster)) %>%
  split(.$trait)

G_GRM <- readRDS(paste0(path_input, "kinships/random_effects/GRMs.rds"))[["Geno"]]

#############################################################
# B) Function for validation
#############################################################

# Helper function (QGG)
run_prediction <- function(data_subset_updated, fixed_model, random_model, GRM_list, validate_list) {
  X_formula <- switch(fixed_model,
    X1 = "~ 1 + G + PC1 + PC2 + PC3",
    X2 = "~ 1 + G + P + PC1 + PC2 + PC3",
    stop("Unsupported fixed_model")
  )

  X_input <- model.matrix(as.formula(X_formula), data = data_subset_updated)

  y <- setNames(data_subset_updated$phenotype, rownames(data_subset_updated))

  model <- tryCatch(
    {
      greml(y = y, X = X_input, GRM = GRM_list, validate = validate_list, ncores = cores)
    },
    error = function(e) {
      cat("Error ", fixed_model, "/", random_model, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(model)) {
    return(NULL)
  }

  temp_output <- tryCatch(
    {
      model$accuracy$Corr %>%
        as.data.frame() %>%
        rownames_to_column("fold") %>%
        mutate("fixed_model" = fixed_model, "random_model" = random_model) %>%
        rename(accuracy = ".")
    },
    error = function(e) {
      cat("Error extracting accuracy for", fixed_model, "/", random_model, ":", conditionMessage(e), "\n")

      return(NULL)
    }
  )
  return(temp_output)
}

# Main function (QGG)
Validation <- function(data_list, G_GRM, validation_scheme) {
  final_output <- data.frame()

  loads <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA.rds"))

  for (current_trait in names(data_list)) {
    print(current_trait)

    meta <- list("trait" = current_trait, "validation_scheme" = validation_scheme, "permuted" = "no", "permutation" = "0")

    data_subset <- data_list[[current_trait]] %>%
      dplyr::mutate(id = as.character(id)) %>%
      drop_na(phenotype, K.Cluster) %>%
      droplevels()

    validate_list <- split(seq_len(nrow(data_subset)), data_subset[[validation_scheme]])

    G_GRM_ext <- G_GRM[data_subset$id, data_subset$id]

    M0 <- list("G" = G_GRM_ext)

    P_GRM_list <- readRDS(paste0(path_input, "kinships/random_effects/GRMs_DA.rds"))

    data_subset_upd <- data_subset %>%
      mutate(.row_order = row_number()) %>%
      left_join(loads %>% distinct(id, G), by = "id") %>%
      arrange(.row_order)

    output <- bind_rows(run_prediction(data_subset_upd, "X1", "M0", M0, validate_list)) %>%
      mutate(!!!meta) %>%
      dplyr::mutate("q_interval" = "baseline")

    final_output <- rbind(final_output, output)

    for (current_interval in names(P_GRM_list)) {
      P_GRM <- P_GRM_list[[current_interval]][["P"]]

      loads_upd <- loads %>% dplyr::filter(q_interval == current_interval)

      data_subset_upd <- data_subset %>%
        dplyr::mutate(.row_order = row_number()) %>%
        left_join(loads_upd, by = "id") %>%
        arrange(.row_order)

      V_list <- list("G" = G_GRM[data_subset_upd$id, data_subset_upd$id], "P" = P_GRM[data_subset_upd$id, data_subset_upd$id])

      output <- bind_rows(
        run_prediction(data_subset_upd, "X1", "M1", V_list, validate_list),
        run_prediction(data_subset_upd, "X2", "M0", M0, validate_list)
      ) %>%
        mutate(!!!meta) %>%
        mutate("q_interval" = current_interval)

      final_output <- rbind(final_output, output)
    }
  }

  saveRDS(final_output, paste0(path_output, "validation/validation_", validation_scheme, "_", submission, ".rds"))
}

#############################################################
# C) Run function
#############################################################

Validation(data_list, G_GRM, "genetic_cluster")

#############################################################
