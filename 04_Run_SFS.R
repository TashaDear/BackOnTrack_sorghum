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
# B) Load input
#############################################################

SFS_mean_list <- readRDS(paste0(path_output, "SFS/SFS_list_input_mean_upd_", submission, ".rds"))

#############################################################
# C) Initilize functions (FastDFE)
#############################################################

BaseInference <- fd$BaseInference
JointInference <- fd$JointInference
Spectra <- fd$Spectra
Spectrum <- fd$Spectrum
shared_params <- fd$SharedParams
Param <- fd$Parametrization
subclass_names <- sapply(reticulate::py_get_attr(Param, "__subclasses__")(), function(cls) reticulate::py_get_attr(cls, "__name__"))

gamma_exp_class <- reticulate::py_get_attr(fd, "GammaExpParametrization")
gamma_exp_instance <- gamma_exp_class()

gamma_dis_class <- reticulate::py_get_attr(fd, "DisplacedGammaParametrization")
gamma_dis_instance <- gamma_dis_class()

input_coefficients <- seq(-100, 100, length.out = 1000)

#############################################################
# D) Marginal SFS function
#############################################################

run_SFS <- function(input, mean_or_random) {
  SFS_parameters_final <- data.frame()
  SFS_spectra_final <- data.frame()
  DFE_discrete_final <- data.frame()

  max_p_b <- 1
  max_S_b <- 100

  for (current_range in names(input)) {
    print(paste0("SNP-category: ", current_range))

    sfs_neutral <- fd$Spectrum(as.integer(as.vector(unlist(input[[current_range]][["neu"]]))))
    sfs_selected <- fd$Spectrum(as.integer(as.vector(unlist(input[[current_range]][["sel"]]))))

    Marginal <- BaseInference(
      sfs_neut = sfs_neutral,
      sfs_sel = sfs_selected,
      n_runs = 100,
      do_bootstrap = TRUE,
      bounds = list(p_b = c(0, max_p_b), S_b = c(0.0001, max_S_b))
    )

    sfs_output <- BaseInference$run(Marginal)

    SFS_param <- Marginal$params_mle %>%
      as.data.frame() %>%
      mutate(
        "analysis" = "Marginal",
        "esm_interval" = current_range,
        "alpha" = Marginal$alpha,
        "theta" = Marginal$theta,
        "llik" = Marginal$likelihood,
        "distribution" = mean_or_random
      )

    SFS_spectra <- data.frame(
      "analysis" = "Marginal",
      "esm_interval" = current_range,
      "bin" = factor(1:20),
      "predicted_sel" = Marginal$sfs_mle$polymorphic,
      "observed_neu" = Marginal$sfs_neut$polymorphic,
      "observed_sel" = Marginal$sfs_sel$polymorphic
    )

    DFE_discrete <- data.frame(
      "analysis" = "Marginal",
      "esm_interval" = current_range,
      "proportion" = BaseInference$get_discretized(Marginal)[[1]],
      "UC" = BaseInference$get_discretized(Marginal)[[2]][1, ],
      "LC" = BaseInference$get_discretized(Marginal)[[2]][2, ],
      "interval" = cut(c(-Inf, -100, -10, -1, 0, 1, Inf)[-1],
        breaks = c(-Inf, -100, -10, -1, 0, 1, Inf), include.lowest = T
      )
    )

    SFS_parameters_final <- rbind(SFS_param, SFS_parameters_final)
    SFS_spectra_final <- rbind(SFS_spectra, SFS_spectra_final)
    DFE_discrete_final <- rbind(DFE_discrete, DFE_discrete_final)
  }

  SFS_DFE_marginal <- list(
    "SFS_param" = SFS_parameters_final,
    "SFS_spectra" = SFS_spectra_final,
    "DFE_dis" = DFE_discrete_final
  )

  saveRDS(SFS_DFE_marginal, paste0(path_output, "SFS/SFS_marginal_analysis_dis_", mean_or_random, "_upd_", submission, ".rds"))
}

#############################################################
# E) Run functions
#############################################################

run_SFS(SFS_mean_list, "mean")

#############################################################
#############################################################
