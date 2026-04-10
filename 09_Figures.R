#############################################################
library("dplyr")
library("data.table")
library("ellipse")
library("emmeans")
library("factoextra")
library("ggpubr")
library("ggsignif")
library("gridExtra")
library("ggExtra")
library("mgcv")
library("multcomp")
library("readr")
library("rlang")
library("viridis")
library("tidyverse")
library("SNPRelate")
library("tibble")
library("tidyr")
library("RColorBrewer")
library("tidygam")
library("tidymv")
library("patchwork")
#############################################################
# A) Paths and input
#############################################################

set.seed(0)

submission <- "submission1"

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/")
path_results <- paste0(path, "results/")

#############################################################
# B) Prepare colorscheme
#############################################################

viridis_colors <- viridis(10)

spectral_colors <- brewer.pal(11, "Spectral")

#############################################################
# C) Functions
#############################################################

sig_levels <- c("permuted", "nonsig.", "sig.", "sig. (Bonferroni)")

cols <- c(
  "[-12.8, -4.8)" = "#440154", "[-4.8, -3.8)" = "#482878", "[-3.8, -3.1)" = "#3E4989", "[-3.1, -2.2)" = "#31688E",
  "[-2.2, -0.9)" = "#26828E", "[-0.9, 0.9)" = "#1F9E89", "[0.9, 2.1)" = "#35B779", "[2.1, 3)" = "#6DCD59",
  "[3, 3.9)" = "#B4DE2C", "[3.9, 11.3)" = "#FDE725", "[3.9, 11.3]" = "#FDE725",
  "neutral" = "grey", "4-fold degenerate" = "grey"
)

trait_recode <- c(
  "BL" = "Terminal Branch Length", "GN" = "Grain Number", "GY" = "Grain Yield", "GW" = "Grain Weight",
  "DTA" = "Days to Anthesis", "PH" = "Plant Height", "PL" = "Panicle Length", "FLH" = "Flag Leaf Height",
  "Starch_adj_for_Protein" = "Starch adjusted", "Protein_adj_for_Starch" = "Protein adjusted"
)
# Function #1
prepare_df <- function(df) {
  df <- df %>% dplyr::mutate(trait = recode(trait, !!!trait_recode))
}

# Function #2
reorder_traits <- function(df) {
  df %>%
    dplyr::mutate(trait = factor(trait, levels = c(
      "Days to Anthesis",
      "Flag Leaf Height", "Panicle Length", "Plant Height", "Terminal Branch Length",
      "Grain Number", "Grain Weight", "Grain Yield",
      "Amylose", "Cal.g", "Fat", "Protein", "Starch"
    ))) %>%
    dplyr::arrange(trait)
}

#############################################################
# D) Load prepared data
#############################################################

Missense_list <- readRDS(paste0(path_input, "/Missense_list.rds"))

Missense_AA <- readRDS(paste0(path_input, "kinships/fixed_effects/Missense_AA_upd.rds")) %>%
  dplyr::filter(Degeneracy != "2") %>%
  dplyr::mutate(afreq_ALT = 1 - afreq, afreq_REF = afreq, afreq_MAF = if_else(afreq_ALT > 0.5, afreq_REF, afreq_ALT)) %>%
  dplyr::select(c(CHROM_POS, AA, Degeneracy, afreq_MAF))

#############################################################

phenotypes <- readRDS(paste0(path_input, "prepared_data/Phenotypes.rds")) %>% do.call(rbind, .)

samples_with_phenotypes <- phenotypes %>%
  distinct(id) %>%
  pull(id)

PCs_geno <- readRDS(paste0(path_input, "prepared_data/PCA_data.rds"))
PCs <- PCs_geno[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::select(c("id", "PC1", "PC2", "PC3"))

passport <- read_csv(file.path(path_input, "SAP/tpj15853-sup-0001-files1.csv")) %>%
  dplyr::rename(id = Taxa) %>%
  dplyr::select(id, K.Cluster, Original_Race) %>%
  dplyr::mutate(BotanicalRace = case_when(
    str_detect(Original_Race, "(?i)verticilliflorum") ~ "Other",
    TRUE ~ Original_Race
  )) %>%
  dplyr::mutate(BotanicalRace = ifelse(is.na(BotanicalRace), "Unknown", BotanicalRace)) %>%
  dplyr::filter(id %in% samples_with_phenotypes, !is.na(K.Cluster)) %>%
  left_join(PCs, by = "id") %>%
  dplyr::mutate(BotanicalRace_upd = str_extract(BotanicalRace, "^[^-]+")) %>%
  dplyr::mutate(BotanicalRace_final = case_when(
    grepl("Bicolor", BotanicalRace, ignore.case = TRUE) ~ "Mixed/Bicolor",
    grepl("-", BotanicalRace) ~ "Mixed", TRUE ~ BotanicalRace
  ))
#############################################################

loads <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA_weighted.rds")) %>% left_join(passport, by = "id")

Scores <- readRDS(paste0(path_input, "kinships/fixed_effects/Score_intervals_final.rds")) %>%
  dplyr::mutate(
    afreq_MAF = if_else(afreq_ALT > 0.5, 1 - afreq_ALT, afreq_ALT),
    afreq_DA = ifelse(AA == ALT, 1 - afreq_ALT, afreq_ALT)
  ) %>%
  dplyr::select(CHROM_POS, AA, ALT, esm_score, esm_score_DA, esm_interval_DA, afreq_DA, afreq_MAF, afreq_ALT)

prop_ALT <- Scores %>%
  group_by(esm_interval_DA) %>%
  summarise(n_total = n(), ALT_is_derived = sum(ALT != AA, na.rm = TRUE), proportion = ALT_is_derived / n_total) %>%
  dplyr::mutate(REF_is_derived = 1 - proportion) %>%
  dplyr::select(esm_interval_DA, ALT_is_derived = proportion, REF_is_derived) %>%
  pivot_longer(cols = c(ALT_is_derived, REF_is_derived), names_to = "Category", values_to = "Proportion")

prop_major <- Scores %>%
  mutate(Category = case_when(
    afreq_DA > 0.5 ~ "Derived is major",
    afreq_DA < 0.5 ~ "Derived is minor", TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Category)) %>%
  group_by(esm_interval_DA, Category) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(Category = factor(Category, levels = c("Derived is major", "Derived is minor")))

#############################################################

LDdecay_ext <- readRDS(paste0(path_input, "LDdecay/LD_prioritized_", submission, ".rds")) %>%
  do.call(rbind, .) %>%
  left_join(Scores, by = "CHROM_POS") %>%
  left_join(Missense_AA, by = "CHROM_POS") %>%
  dplyr::filter(adjustment == "Q+K") %>%
  dplyr::mutate(
    SNP_category = if_else(is.na(esm_interval_DA), "neutral", esm_interval_DA),
    degeneracy = paste0("fold", Degeneracy)
  )

#############################################################
# E) Load model output
#############################################################

Mean_partition <- rbind(
  readRDS(paste0(path_output, "mean_partition/Mean_partition_", submission, ".rds")),
  readRDS(paste0(path_output, "mean_partition/Mean_partition_permuted_upd_", submission, ".rds"))
) %>% prepare_df()

Variance_components <- readRDS(paste0(path_output, "variance_partition/Variance_components_", submission, ".rds")) %>% prepare_df()
Variance_partition <- rbind(
  readRDS(paste0(path_output, "variance_partition/Likelihoods_", submission, ".rds")),
  readRDS(paste0(path_output, "variance_partition/Likelihoods_permuted_", submission, ".rds"))
) %>% prepare_df()

Validation <- rbind(
  readRDS(paste0(path_output, "validation/validation_genetic_cluster_", submission, ".rds")),
  readRDS(paste0(path_output, "validation/validation_genetic_cluster_permuted_", submission, ".rds"))
) %>% prepare_df()

SFS_DFE_marginal <- readRDS(paste0(path_output, "SFS/SFS_marginal_analysis_", submission, "_new_version.rds"))

SFS_DFE_joint <- readRDS(paste0(path_output, "SFS/SFS_joint_analysis_", submission, ".rds"))

#############################################################
# F) Distribution of REF/ALT as derived alleles in different SNP-categories
#############################################################

fill_colors <- c(
  "ALT_is_derived" = spectral_colors[5],
  "REF_is_derived" = spectral_colors[9],
  "Derived is major" = spectral_colors[4],
  "Derived is minor" = spectral_colors[7],
  "score (ALT allele)" = spectral_colors[5],
  "score (Derived allele)" = spectral_colors[8],
  "score (REF allele)" = spectral_colors[9]
)

Bar_plot <- ggplot(prop_ALT, aes(x = esm_interval_DA, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.75) +
  scale_fill_manual(values = fill_colors, labels = c("ALT is derived", "REF is derived")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Evolutionary score interval", y = "Proportion", fill = "") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.justification = "center",
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 5),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

Bar_v2_plot <- ggplot(prop_major, aes(x = esm_interval_DA, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.75) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Evolutionary score interval", y = "Proportion", fill = "") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.justification = "center",
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 5),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

Scores_upd <- Scores %>%
  dplyr::mutate(esm_score_ALT = esm_score, esm_score_REF = -esm_score) %>%
  pivot_longer(cols = c(esm_score_ALT, esm_score_REF, esm_score_DA), names_to = "Polarization", values_to = "Score") %>%
  dplyr::mutate(
    Polarization = recode(Polarization,
      "esm_score_ALT" = "score (ALT allele)",
      "esm_score_DA" = "score (Derived allele)",
      "esm_score_REF" = "score (REF allele)"
    ),
    Polarization = factor(Polarization,
      levels = c("score (ALT allele)", "score (REF allele)", "score (Derived allele)")
    )
  )

Histogram_pol_plot <- ggplot(Scores_upd, aes(x = Score, fill = Polarization)) +
  geom_histogram(alpha = 1, bins = 100, position = "identity", color = NA) +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.7, linetype = "dashed") +
  scale_fill_manual(values = fill_colors) +
  labs(x = "Evolutionary scores", y = "") +
  theme_bw() +
  theme(
    legend.position = c(0.01, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.margin = margin(t = 5, r = 15, b = 5, l = 5),
    legend.direction = "vertical"
  ) +
  xlim(-20, 20)

png(paste0(path_results, "figures/evolutionary_scores/Histogram_scores_polarization_updated_", submission, ".png"),
  width = 7, height = 3.5, units = "in", res = current_dpi
)
Histogram_pol_plot_upd <- Histogram_pol_plot + annotation_custom(ggplotGrob(Bar_v2_plot), xmin = 5, xmax = 22, ymin = 500, ymax = 1900)
print(Histogram_pol_plot_upd)
dev.off()

#############################################################
# G) Load of prioritized alleles
#############################################################

loads <- loads %>% dplyr::mutate(K.Cluster = as.factor(K.Cluster), BotanicalRace_final = as.factor(BotanicalRace_final))

overall_mean <- mean(loads$P, na.rm = TRUE)

loads_v1 <- loads %>% dplyr::filter(!is.na(K.Cluster))
load_mean_v1 <- mean(loads_v1$P, na.rm = TRUE)

loads_v2 <- loads %>% dplyr::filter(!BotanicalRace_final %in% c("Mixed", "Other", "Unknown"), !is.na(BotanicalRace_final))
load_mean_v2 <- mean(loads_v2$P, na.rm = TRUE)

#############################################################

# Loads- genetic cluster
anova_model_v1 <- aov(P ~ K.Cluster, data = loads_v1)
summary(anova_model_v1)

library(emmeans)
emmeans(anova_model_v1, pairwise ~ K.Cluster)


tukey_results_v1 <- TukeyHSD(anova_model_v1)$K.Cluster %>%
  as.data.frame() %>%
  rownames_to_column("comparison") %>%
  dplyr::mutate(across(-comparison, ~ round(.x, 5)))

table_plot_v1 <- tableGrob(tukey_results_v1, rows = NULL)
ggsave(paste0(path_results, "figures/evolutionary_scores/Genetic_cluster_load.pdf"), plot = table_plot_v1, width = 8, height = 6)

cluster_counts_v1 <- loads_v1 %>%
  group_by(K.Cluster) %>%
  summarise(
    n = n(), y_pos = max(P, na.rm = TRUE) + 0.05 * abs(max(P, na.rm = TRUE)),
    label = paste0("n = ", n), .groups = "drop"
  )

Allele_load_v1_plot <- ggplot(loads_v1, aes(x = K.Cluster, y = P, fill = K.Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.05, color = "black", size = 0.5) +
  geom_hline(yintercept = load_mean_v1, linetype = "dashed", color = "black") +
  scale_fill_viridis_d(name = "Genetic cluster") +
  geom_text(data = cluster_counts_v1, aes(x = as.factor(K.Cluster), y = 3500, label = label), inherit.aes = FALSE, size = 3) +
  labs(x = "", y = "Weighted allele load") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    strip.text = element_text(size = 12),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  guides(fill = guide_legend(title = "Genetic Cluster", nrow = 1)) +
  ylim(0, 4000)

# Loads- Botanical Race
all_races <- levels(factor(loads_v2$BotanicalRace_final))

color_vec <- setNames(
  c(
    brewer.pal(length(all_races), "Spectral"),
    rep("grey70", length(setdiff(unique(passport$BotanicalRace_final), all_races)))
  ),
  c(all_races, setdiff(unique(passport$BotanicalRace_final), all_races))
)

anova_model_v2 <- aov(P ~ BotanicalRace_final, data = loads_v2)
tukey_results_v2 <- TukeyHSD(anova_model_v2)$BotanicalRace_final %>%
  as.data.frame() %>%
  rownames_to_column("comparison") %>%
  dplyr::mutate(across(-comparison, ~ round(.x, 5)))

table_plot_v2 <- tableGrob(tukey_results_v2, rows = NULL)
ggsave(paste0(path_results, "figures/evolutionary_scores/Botanical_Race_load.pdf"), plot = table_plot_v2, width = 10, height = 4)

cluster_counts_v2 <- loads_v2 %>%
  group_by(BotanicalRace_final) %>%
  summarise(
    n = n(), y_pos = max(P, na.rm = TRUE) + 0.05 * abs(max(P, na.rm = TRUE)),
    label = paste0("n = ", n), .groups = "drop"
  )

Allele_load_v2_plot <- ggplot(loads_v2, aes(x = as.factor(BotanicalRace_final), y = P, fill = as.factor(BotanicalRace_final))) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.05, color = "black", size = 0.5) +
  geom_hline(yintercept = load_mean_v2, linetype = "dashed", color = "black") +
  scale_fill_manual(values = color_vec, name = "Botanical Race") +
  geom_text(
    data = cluster_counts_v2, aes(x = as.factor(BotanicalRace_final), y = 3500, label = label),
    inherit.aes = FALSE, size = 3
  ) +
  labs(x = "", y = "Weighted allele load") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    strip.text = element_text(size = 12),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  guides(fill = guide_legend(title = "Botanical Race", nrow = 1)) +
  ylim(0, 4000)

#############################################################
# H) PCA plots
#############################################################

# PCA- genetic cluster
PCA_cluster_plot <- passport %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(K.Cluster)), size = 1.5, alpha = 1) +
  labs(
    color = "Genetic cluster       ",
    x = sprintf("PC1 (%.2f%%)", (PCs_geno$sdev[1]^2 / sum(PCs_geno$sdev^2) * 100)),
    y = sprintf("PC2 (%.2f%%)", (PCs_geno$sdev[2]^2 / sum(PCs_geno$sdev^2) * 100))
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.05, 0.025),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.5),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  scale_color_viridis_d(option = "D") +
  ylim(-2000, 2000) +
  xlim(-1500, 1500) +
  guides(color = guide_legend(nrow = 2))
ggsave(paste0(path_results, "figures/PCA/PCA_genetic_cluster_", submission, ".png"),
  plot = PCA_cluster_plot, width = 10, height = 5, dpi = current_dpi
)

# PCA- Botanical Race
PCA_race_plot <- passport %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(BotanicalRace_final)), size = 1.5, alpha = 1) +
  labs(
    color = "Botanical Race       ",
    x = sprintf("PC1 (%.2f%%)", (PCs_geno$sdev[1]^2 / sum(PCs_geno$sdev^2) * 100)),
    y = sprintf("PC2 (%.2f%%)", (PCs_geno$sdev[2]^2 / sum(PCs_geno$sdev^2) * 100))
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.05, 0.025),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.5),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0.1, "cm")
  ) +
  scale_color_manual(values = color_vec, name = "Botanical Race") +
  ylim(-2000, 2000) +
  xlim(-1500, 1500) +
  guides(color = guide_legend(nrow = 2))
ggsave(paste0(path_results, "figures/PCA/PCA_botanical_race_", submission, ".png"),
  plot = PCA_race_plot, width = 10, height = 5, dpi = current_dpi
)

#############################################################

png(paste0(path_results, "figures/PCA/PCA_genetic_cluster_extended_", submission, ".png"), width = 7, height = 5, units = "in", res = current_dpi)
PCA_cluster_plot_upd <- PCA_cluster_plot + annotation_custom(ggplotGrob(Allele_load_v1_plot), xmin = -1660, xmax = 450, ymin = 500, ymax = 2140)
print(PCA_cluster_plot_upd)
dev.off()

png(paste0(path_results, "figures/PCA/PCA_botanical_race_extended_", submission, ".png"), width = 7, height = 5, units = "in", res = current_dpi)
PCA_race_plot_upd <- PCA_race_plot + annotation_custom(ggplotGrob(Allele_load_v2_plot), xmin = -1660, xmax = 450, ymin = 500, ymax = 2140)
print(PCA_race_plot_upd)
dev.off()

#############################################################
# I) Phenotypes plot
#############################################################

Phenotypes_df <- phenotypes %>%
  dplyr::mutate(trait = recode(trait,
    "BL" = "Terminal Branch Length", "GN" = "Grain Number", "GY" = "Grain Yield",
    "GW" = "Grain Weight", "DTA" = "Days to Anthesis", "PH" = "Plant Height",
    "PL" = "Panicle Length", "FLH" = "Flag Leaf Height"
  )) %>%
  dplyr::filter(!trait %in% c("Plant Height"), !is.na(phenotype)) %>%
  reorder_traits()

Phenotypes_plot <- ggplot(Phenotypes_df, aes(x = phenotype)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black", alpha = 0.6) +
  facet_wrap(~trait, ncol = 6, nrow = 2, scales = "free_x") +
  theme_bw() +
  labs(y = "Frequency", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#############################################################
# J) Combine plots
#############################################################

S1_AB <- ggarrange(PCA_race_plot_upd, Histogram_pol_plot_upd, ncol = 2, labels = c("A", "B"), widths = c(1, 1.25))
S1_C <- ggarrange(Phenotypes_plot, ncol = 1, labels = c("C"))
S1_combined <- ggarrange(S1_AB, S1_C, nrow = 2, heights = c(1, 1.10))

#############################################################
# K) Mean partition
#############################################################

color_to_use <- viridis_colors[3]

Mean_partition_upd <- Mean_partition %>%
  dplyr::mutate(zscore = effect / se)

Mean_part_pempirical <- Mean_partition_upd %>%
  group_by(trait, q_interval) %>%
  summarise(
    observed_zscore = zscore[permuted == "no"][1],
    r = {
      obs <- zscore[permuted == "no"][1]
      permuted_yes <- zscore[permuted == "yes"]
      if (obs >= 0) {
        sum(permuted_yes >= obs)
      } else {
        sum(permuted_yes <= obs)
      }
    },
    n = sum(permuted == "yes"),
    p_empirical = (r + 1) / (n + 1), .groups = "drop"
  ) %>%
  ungroup()

Mean_part_combined <- Mean_partition_upd %>%
  left_join(Mean_part_pempirical, by = c("trait", "q_interval")) %>%
  dplyr::mutate(
    significance =
      case_when(
        permuted == "yes" ~ "permuted",
        permuted == "no" & p_empirical <= (0.05 / 10) & sig <= 0.05 ~ "sig. (Bonferroni)",
        permuted == "no" & p_empirical <= 0.05 & p_empirical > (0.05 / 10) & sig <= 0.05 ~ "sig.",
        permuted == "no" ~ "nonsig."
      ),
    significance = factor(significance, levels = c("permuted", "nonsig.", "sig.", "sig. (Bonferroni)"))
  ) %>%
  dplyr::mutate(
    esm_start = as.numeric(str_extract(esm_interval, "(?<=\\[)[^,]+")),
    esm_interval = factor(esm_interval, levels = unique(esm_interval[order(esm_start)])),
    esm_interval = recode(esm_interval, "[3.9, 11.3)" = "[3.9, 11.3]")
  ) %>%
  dplyr::filter(trait != "Plant Height") %>%
  dplyr::mutate(analysis = "Mean partition") %>%
  reorder_traits()

Mean_binned_plot <- ggplot(Mean_part_combined, aes(x = esm_interval, y = zscore)) +
  geom_point(
    data = filter(Mean_part_combined, significance == "permuted"), aes(shape = significance),
    color = "grey60", fill = NA, size = 1, stroke = 0.5, alpha = 0.5
  ) +
  geom_point(
    data = filter(Mean_part_combined, significance == "nonsig."), aes(shape = significance),
    color = color_to_use, fill = NA, size = 1.3, stroke = 1.2
  ) +
  geom_point(
    data = filter(Mean_part_combined, significance == "sig."), aes(shape = significance),
    color = color_to_use, fill = color_to_use, size = 4, stroke = 0.2
  ) +
  geom_point(
    data = filter(Mean_part_combined, significance == "sig. (Bonferroni)"), aes(shape = significance),
    color = color_to_use, fill = color_to_use, size = 4, stroke = 0.2
  ) +
  geom_line(
    data = filter(Mean_part_combined, permuted == "no"), aes(group = trait),
    color = "black", linewidth = 0.4, alpha = 0.6
  ) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_grid(trait ~ analysis) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.box = "vertical",
    legend.spacing.x = unit(6, "mm"),
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.key = element_rect(fill = "transparent", color = NA)
  ) +
  scale_shape_manual(
    name = "Empirical p-value ",
    values = c("permuted" = 1, "nonsig." = 21, "sig." = 21, "sig. (Bonferroni)" = 23),
    labels = c("Permuted", "P-value > 0.05", "P-value < 0.05", "P-value < 0.005")
  ) +
  scale_fill_manual(
    name = "Empirical p-value ",
    values = c("permuted" = NA, "nonsig." = NA, "sig." = color_to_use, "sig. (Bonferroni)" = color_to_use)
  ) +
  scale_color_manual(
    name = "Empirical p-value ",
    values = c("permuted" = "grey60", "nonsig." = color_to_use, "sig." = color_to_use)
  ) +
  guides(shape = guide_legend(
    override.aes =
      list(
        size = c(2, 2, 4, 5), fill = c(NA, NA, "grey", "grey"), color = c("grey60", "grey", "grey", "grey"),
        stroke = c(0.5, 1.2, 0.2, 0.2)
      )
  ), fill = "none", color = "none", size = "none") +
  labs(y = "Z-scores", x = "") +
  scale_x_discrete(drop = TRUE)

# Table
Mean_part_table <- Mean_part_combined %>%
  dplyr::filter(permuted == "no", significance != "nonsig.") %>%
  dplyr::rename(p_wald = sig) %>%
  dplyr::select(c("trait", "esm_interval", "effect", "se", "p_wald", "p_empirical")) %>%
  mutate(across(where(is.numeric), ~ signif(., 5))) %>%
  dplyr::mutate("analysis" = "Mean partition") %>%
  tableGrob(rows = NULL)

#############################################################
# L) Variance partition- components
#############################################################

color_vec <- c(
  "G" = viridis_colors[4],
  "P" = viridis_colors[9]
)

Scores_DA_distinct <- Scores %>%
  distinct(esm_interval_DA) %>%
  arrange(esm_interval_DA) %>%
  dplyr::mutate(q_interval = paste0("[", seq(0.0, 0.9, 0.1), ", ", seq(0.1, 1, 0.1), ")"))

Variance_components_upd <- Variance_components %>%
  dplyr::filter(model == "M1", component != "error") %>%
  reorder_traits() %>%
  left_join(Scores_DA_distinct, by = "q_interval") %>%
  mutate(analysis = "Variance components") %>%
  dplyr::filter(trait != "Plant Height")

variance_plot <- ggplot(Variance_components_upd, aes(x = esm_interval_DA, y = sigma, color = component)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(group = component)) +
  facet_wrap(~trait, ncol = 4, nrow = 3, scales = "free_y") +
  scale_color_manual(values = color_vec) +
  theme_bw() +
  labs(x = "", y = "Variance", color = "Component") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(2, "lines"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

#############################################################

cols <- c(
  "[-12.8, -4.8)" = "#440154", "[-4.8, -3.8)" = "#482878", "[-3.8, -3.1)" = "#3E4989", "[-3.1, -2.2)" = "#31688E", "[-2.2, -0.9)" = "#26828E",
  "[-0.9, 0.9)" = "#1F9E89", "[0.9, 2.1)" = "#35B779", "[2.1, 3)" = "#6DCD59", "[3, 3.9)" = "#B4DE2C", "[3.9, 11.3)" = "#FDE725",
  "[3.9, 11.3]" = "#FDE725", "G" = "grey"
)

G_scaling <- readRDS(paste0(path_input, "kinships/random_effects/GRMs.rds"))[["Scaling"]] %>%
  dplyr::filter(GRM == "Geno") %>%
  pull(Scaling)

GRM_P_list <- readRDS(paste0(path_input, "kinships/random_effects/GRMs_DA.rds"))

scaling_df <- lapply(names(GRM_P_list), function(q) {
  GRM_P_list[[q]]$Scaling %>% dplyr::mutate(q_interval = q)
}) %>%
  bind_rows() %>%
  rename(P_scaling = Scaling)

Variance_components_sub <- Variance_components %>%
  dplyr::filter(model == "M1") %>%
  left_join(scaling_df, by = "q_interval") %>%
  mutate(sigma_rescaled = case_when(
    component == "G" ~ sigma / G_scaling,
    component == "P" ~ sigma / P_scaling, TRUE ~ NA_real_
  )) %>%
  dplyr::filter(component != "error") %>%
  left_join(Scores_DA_distinct, by = "q_interval") %>%
  group_by(trait, q_interval) %>%
  mutate(prop_var = sigma / sum(sigma, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-c("GRM", "model", "P_scaling", "permuted", "permutation")) %>%
  mutate(esm_interval_upd = if_else(component == "G", "G", esm_interval_DA))

SXX_plot <- ggplot(Variance_components_sub, aes(x = trait, y = prop_var, fill = esm_interval_upd)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "Proportion of genetic variance") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~esm_interval_DA, nrow = 2, ncol = 5, scales = "free_x")

#############################################################
# M) Variance partition- LLR
#############################################################

color_to_use <- viridis_colors[7]

M0_df <- Variance_partition %>%
  dplyr::filter(model == "M0", permuted == "no") %>%
  dplyr::select(c("trait", "llik")) %>%
  dplyr::rename("M0_llik" = "llik")

LLRs_combined <- Variance_partition %>%
  left_join(M0_df, by = "trait") %>%
  dplyr::mutate(
    LLR_diff = 2 * (llik - M0_llik),
    pval = if_else(permuted == "no", round(pchisq(LLR_diff, df = 1, lower.tail = FALSE), 4), NA_real_),
    permutation = if_else(permuted == "no", "none", permutation)
  ) %>%
  dplyr::filter(permuted == "yes" | (permuted == "no" & q_interval != "baseline")) %>%
  dplyr::select(-c(M0_llik, llik)) %>%
  mutate(analysis = "Variance partition", sig = pval)

LLRs_pempirical <- LLRs_combined %>%
  group_by(trait, q_interval) %>%
  summarise(
    observed_LLR = LLR_diff[permuted == "no"][1],
    r = sum(LLR_diff[permuted == "yes"] >= observed_LLR),
    n = sum(permuted == "yes"),
    p_empirical = (r + 1) / (n + 1), .groups = "drop"
  )

Var_part_combined <- LLRs_combined %>%
  left_join(LLRs_pempirical, by = c("trait", "q_interval")) %>%
  dplyr::mutate(
    significance = case_when(
      permuted == "yes" ~ "permuted",
      permuted == "no" &
        p_empirical <= (0.05 / 10) & sig <= 0.05 ~ "sig. (Bonferroni)",
      permuted == "no" &
        p_empirical <= 0.05 & p_empirical > (0.05 / 10) & sig <= 0.05 ~ "sig.",
      permuted == "no" ~ "nonsig."
    ),
    significance = factor(significance, levels = c("permuted", "nonsig.", "sig.", "sig. (Bonferroni)")),
    point_size = ifelse(permuted == "yes", 0.5, 5)
  ) %>%
  dplyr::mutate(significance = factor(significance, levels = c("permuted", "nonsig.", "sig.", "sig. (Bonferroni)"))) %>%
  left_join(Scores_DA_distinct, by = "q_interval") %>%
  dplyr::mutate(esm_interval = recode(esm_interval_DA, "[3.9, 11.3)" = "[3.9, 11.3]")) %>%
  dplyr::filter(trait != "Plant Height") %>%
  dplyr::mutate(analysis = "Variance partition") %>%
  reorder_traits()

Var_binned_plot <- ggplot(Var_part_combined, aes(x = esm_interval, y = LLR_diff)) +
  geom_point(
    data = filter(Var_part_combined, significance == "permuted"), aes(shape = significance),
    color = "grey60", fill = NA, size = 1, stroke = 0.5, alpha = 0.5
  ) +
  geom_line(
    data = filter(Var_part_combined, permuted == "no"), aes(group = trait),
    color = "black", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = filter(Var_part_combined, significance == "nonsig."), aes(shape = significance),
    color = color_to_use, fill = NA, size = 1.3, stroke = 1.2
  ) +
  geom_point(
    data = filter(Var_part_combined, significance == "sig."), aes(shape = significance),
    color = color_to_use, fill = color_to_use, size = 4, stroke = 0.2
  ) +
  geom_point(
    data = filter(Var_part_combined, significance == "sig. (Bonferroni)"), aes(shape = significance),
    color = color_to_use, fill = color_to_use, size = 4, stroke = 0.2
  ) +
  geom_hline(yintercept = qchisq(0.95, df = 1), linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_grid(trait ~ analysis) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box = "vertical",
    legend.spacing.x = unit(6, "mm"),
    legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_shape_manual(
    name = "Empirical p-value ",
    values = c("permuted" = 1, "nonsig." = 21, "sig." = 21, "sig. (Bonferroni)" = 23),
    labels = c("Permuted", "P-value > 0.05", "P-value < 0.05", "P-value < 0.005")
  ) +
  scale_fill_manual(
    name = "Empirical p-value ",
    values = c("permuted" = NA, "nonsig." = NA, "sig." = color_to_use, "sig. (Bonferroni)" = color_to_use)
  ) +
  scale_color_manual(
    name = "Empirical p-value ",
    values = c("permuted" = "grey60", "nonsig." = color_to_use, "sig." = color_to_use)
  ) +
  guides(shape = guide_legend(
    override.aes =
      list(
        size = c(2, 2, 4, 5), fill = c(NA, NA, "grey", "sig. (Bonferroni)" = "grey"),
        color = c("grey60", "grey", "grey", "grey"),
        stroke = c(0.5, 1.2, 0.2, 0.2)
      )
  ), fill = "none", color = "none", size = "none") +
  labs(y = "LLRs", x = "") +
  scale_x_discrete(drop = TRUE)

# Table
Var_part_table <- Var_part_combined %>%
  dplyr::filter(permuted == "no", significance != "nonsig.") %>%
  dplyr::rename(p_analytical = pval) %>%
  dplyr::select(c("trait", "esm_interval", "LLR_diff", "p_analytical", "p_empirical")) %>%
  dplyr::mutate("analysis" = "Variance partition") %>%
  tableGrob(rows = NULL)

#############################################################
# N) Validation
#############################################################

color_to_use_v2 <- c(
  "Mean partition" = viridis_colors[3],
  "Variance partition" = viridis_colors[7]
)

Validation_upd <- Validation %>%
  group_by(trait, random_model, fixed_model, validation_scheme, permutation, permuted, q_interval) %>%
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    se_accuracy = sd(accuracy, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup() %>%
  dplyr::filter(validation_scheme == "genetic_cluster")

M0_PA <- Validation_upd %>%
  dplyr::filter(random_model == "M0" & fixed_model == "X1", permuted == "no") %>%
  dplyr::select(c("trait", "mean_accuracy", "se_accuracy")) %>%
  dplyr::rename("baseline" = "mean_accuracy", "baseline (se)" = "se_accuracy")

Validation_final <- Validation_upd %>%
  left_join(M0_PA, by = c("trait")) %>%
  dplyr::filter(!(random_model == "M0" & fixed_model == "X1")) %>%
  dplyr::mutate(
    improv = mean_accuracy - baseline,
    partition = case_when(
      fixed_model == "X2" ~ "Mean partition",
      fixed_model != "X2" ~ "Variance partition", TRUE ~ NA_character_
    )
  )

Validation_pempirical <- Validation_final %>%
  group_by(trait, partition, q_interval) %>%
  summarise(
    observed_improv = improv[permuted == "no"][1],
    r = sum(improv[permuted == "yes"] >= observed_improv),
    n = sum(permuted == "yes"),
    p_empirical = (r + 1) / (n + 1), .groups = "drop"
  )

Validation_combined <- Validation_final %>%
  left_join(Validation_pempirical,
    by = c("trait", "partition", "q_interval")
  ) %>%
  dplyr::mutate(
    significance =
      case_when(
        permuted == "yes" ~ "permuted",
        permuted == "no" & p_empirical <= 0.05 & improv > 0 ~ "sig.",
        permuted == "no" & p_empirical <= (0.05 / 10) & improv > 0 ~ "sig. (Bonferroni)",
        permuted == "no" & p_empirical > 0.05 | improv <= 0 ~ "nonsig."
      ),
    point_size = ifelse(permuted == "yes", 0.5, 5),
    analysis = "Validation"
  ) %>%
  left_join(Scores_DA_distinct, by = "q_interval") %>%
  rename("esm_interval" = "esm_interval_DA") %>%
  dplyr::mutate(
    esm_start = as.numeric(str_extract(esm_interval, "(?<=\\[)[^,]+")),
    esm_interval = factor(esm_interval, levels = unique(esm_interval[order(esm_start)])),
    esm_interval = recode(esm_interval, "[3.9, 11.3)" = "[3.9, 11.3]")
  ) %>%
  dplyr::filter(trait != "Plant Height")

missing_levels <- setdiff(sig_levels, unique(Validation_combined$significance))
if (length(missing_levels) > 0) {
  dummy_rows <- tibble::tibble(
    esm_interval = unique(Validation_combined$esm_interval)[1],
    improv = NA,
    trait = unique(Validation_combined$trait)[1],
    permuted = "no",
    significance = factor(missing_levels, levels = sig_levels)
  )

  Validation_combined <- dplyr::bind_rows(Validation_combined, dummy_rows)
}

Validation_combined <- Validation_combined %>%
  dplyr::mutate(group = paste0(trait, partition), analysis = "Validation") %>%
  dplyr::filter(permuted != "yes")

Val_binned_plot <- ggplot(Validation_combined, aes(x = esm_interval, y = improv)) +
  geom_line(
    data = filter(Validation_combined, permuted == "no"), aes(group = partition, color = partition),
    linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = filter(Validation_combined, significance == "nonsig."),
    aes(shape = significance, color = partition), fill = NA, size = 1.5, stroke = 1.2
  ) +
  geom_point(
    data = filter(Validation_combined, significance %in% c("sig.", "sig. (Bonferroni)")),
    aes(shape = significance, color = partition, fill = partition, size = significance), stroke = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_grid(trait ~ analysis) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.justification = "center",
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  scale_color_manual(name = "Partition", values = color_to_use_v2, guide = guide_legend(order = 1), na.translate = FALSE) +
  scale_fill_manual(name = "Partition", values = color_to_use_v2, guide = guide_legend(order = 1), na.translate = FALSE) +
  scale_shape_manual(
    name = "Empirical p-value  ",
    values = c("nonsig." = 21, "sig." = 21, "sig. (Bonferroni)" = 23), labels = c("P-value > 0.05", "P-value < 0.05", "P-value < 0.005"),
    guide = guide_legend(order = 2, override.aes = list(
      size = c(3, 4, 4),
      fill = c(alpha("#4D4D4D", 0.5), "#4D4D4D", "#4D4D4D"),
      color = c("#1A1A1A", "#4D4D4D", "#4D4D4D"), stroke = c(2, 0.2, 0.2)
    ))
  ) +
  scale_size_manual(values = c("permuted" = 1, "nonsig." = 1.5, "sig." = 4, "sig. (Bonferroni)" = 4), guide = "none") +
  labs(y = "Improvement in PA", x = "") +
  scale_x_discrete(drop = TRUE) +
  ylim(-0.15, 0.15) +
  theme(
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0)
  )

# Table
Val_part_table <- Validation_combined %>%
  dplyr::filter(permuted == "no", !is.na(improv), significance != "nonsig.") %>%
  dplyr::select(c(
    "trait", "esm_interval", "partition",
    "baseline", "mean_accuracy", "se_accuracy", "improv", "p_empirical", "analysis"
  )) %>%
  tableGrob(rows = NULL)

#############################################################
# O) Combined plots
#############################################################

# Data list
Combined_data_list <- list("Mean partition" = Mean_part_combined, "Variance partition" = Var_part_combined, "Validation" = Validation_combined)
saveRDS(Combined_data_list, paste0(path_output, "Combined_data_", submission, ".rds"))

# Combined plots
Mean_part_plot_upd <- Mean_binned_plot + theme(legend.position = "bottom") +
  theme(strip.text.x = element_text(face = "bold", size = 11))

Var_part_plot_upd <- Var_binned_plot + theme(legend.position = "none") +
  theme(strip.text.x = element_text(face = "bold", size = 11))

Val_binned_plot <- Val_binned_plot + theme(legend.position = "none") +
  theme(strip.text.x = element_text(face = "bold", size = 11))

plots <- ggarrange(Mean_part_plot_upd, Var_part_plot_upd, Val_binned_plot, ncol = 3, labels = c("A", "B", "C"), common.legend = T, legend = "bottom")

#############################################################
# P) LDdecay: Prepare input
#############################################################

SNPs_after_filtering <- Missense_AA %>%
  dplyr::filter(afreq_MAF > 0.05) %>%
  pull(CHROM_POS)

LDdecay_ext_final <- LDdecay_ext %>%
  mutate(SNP_category = factor(recode(SNP_category, "[3.9, 11.3)" = "[3.9, 11.3]", "neutral" = "4-fold degenerate"),
    levels = c(
      "4-fold degenerate", "[-12.8, -4.8)", "[-4.8, -3.8)", "[-3.8, -3.1)", "[-3.1, -2.2)",
      "[-2.2, -0.9)", "[-0.9, 0.9)", "[0.9, 2.1)", "[2.1, 3)", "[3, 3.9)", "[3.9, 11.3]"
    )
  )) %>%
  dplyr::filter(CHROM_POS %in% SNPs_after_filtering, !is.na(SNP_category))

#############################################################
# Q) LDdecay: Prepare models
#############################################################

# Determine best model
LDdecay_model_log <- glm(r2 ~ D + D:SNP_category, data = LDdecay_ext_final, family = Gamma(link = "log"))
LDdecay_model_inverse <- glm(r2 ~ D + D:SNP_category, data = LDdecay_ext_final, family = Gamma(link = "inverse"))
AIC(LDdecay_model_log, LDdecay_model_inverse)

# Determine if interaction-term is significant
LDdecay_model_inverse_A <- glm(r2 ~ D, data = LDdecay_ext_final, family = Gamma(link = "inverse"))
LDdecay_model_inverse_B <- glm(r2 ~ D + SNP_category, data = LDdecay_ext_final, family = Gamma(link = "inverse"))
LDdecay_model_inverse_C <- glm(r2 ~ D + SNP_category + D:SNP_category, data = LDdecay_ext_final, family = Gamma(link = "inverse"))

intercept_table <- anova(LDdecay_model_inverse_A, LDdecay_model_inverse_B, test = "Chisq") %>%
  as.data.frame() %>%
  mutate("analysis" = "LD-intercept")

slope_table <- anova(LDdecay_model_inverse_B, LDdecay_model_inverse_C, test = "Chisq") %>%
  as.data.frame() %>%
  mutate("analysis" = "slopes")

LD_anova <- rbind(intercept_table, slope_table) %>% tableGrob(rows = NULL)

#############################################################
# Q) LDdecay: Prepare plot
#############################################################

current_model <- LDdecay_model_inverse_C

pred_data <- LDdecay_ext_final %>%
  distinct(SNP_category) %>%
  expand_grid(D = seq(0, max(LDdecay_ext_final$D), length.out = 500)) %>%
  dplyr::mutate(
    r2_pred = predict(current_model, newdata = ., type = "response"),
    category = factor(SNP_category, levels = levels(LDdecay_ext_final$category))
  )

linetype_vals <- setNames(
  ifelse(levels(LDdecay_ext_final$SNP_category) == "4-fold degenerate", "dashed", "solid"),
  levels(LDdecay_ext_final$SNP_category)
)

pred_data <- pred_data %>% arrange(SNP_category, D)

LD_decay_plot <- ggplot(LDdecay_ext_final, aes(x = D, y = r2, color = SNP_category)) +
  geom_point(alpha = 0) +
  geom_line(data = pred_data, aes(x = D, y = r2_pred, color = SNP_category, linetype = SNP_category), linewidth = 0.5) +
  geom_hline(yintercept = 0.2, color = "red", linetype = "dashed") +
  scale_color_manual(values = cols) +
  labs(x = "Distance (bp)", y = expression(LD ~ (r^2)), color = "SNP category") +
  theme_bw() +
  theme(legend.position = "left") +
  xlim(0, 50000) +
  ylim(0, 0.4) +
  scale_linetype_manual(values = linetype_vals, guide = "none") +
  guides(color = guide_legend(nrow = 11))

ggsave(paste0(path_results, "figures/LDdecay/LDdecay_", submission, ".png"), plot = LD_decay_plot, width = 8, height = 4, dpi = 500)

#############################################################
# S) SFS analysis
#############################################################

current_path <- paste0(path_results, "figures/SFS/")

bin_order <- c(
  "[-12.8, 11.3]", "[-12.8, -4.8)", "[-4.8, -3.8)", "[-3.8, -3.1)", "[-3.1, -2.2)",
  "[-2.2, -0.9)", "[-0.9, 0.9)", "[0.9, 2.1)", "[2.1, 3)", "[3, 3.9)", "[3.9, 11.3)", "[3.9, 11.3]"
)

prepare_SFS <- function(df) {
  df <- df %>%
    dplyr::mutate(
      esm_interval = recode(esm_interval, "Selected" = "[-12.8, 11.3]"),
      esm_interval = factor(esm_interval, levels = bin_order)
    )
}

SFS <- SFS_DFE_marginal[["SFS_spectra"]] %>%
  pivot_longer(cols = c(predicted_sel, observed_neu, observed_sel), names_to = "Spectra", values_to = "Count") %>%
  group_by(esm_interval, Spectra) %>%
  dplyr::mutate(Proportion = Count / sum(Count), interval = bin) %>%
  ungroup() %>%
  prepare_SFS()

SFS_sel <- SFS %>% dplyr::filter(Spectra == "predicted_sel")

DFE_discrete <- SFS_DFE_marginal[["DFE_dis"]] %>% prepare_SFS()

SFS_param <- SFS_DFE_marginal[["SFS_param"]] %>%
  dplyr::mutate(
    lower = as.numeric(str_extract(esm_interval, "-?\\d+\\.?\\d*(?=,)")),
    order_group = ifelse(esm_interval == "Selected", 0, 1)
  ) %>%
  arrange(order_group, lower) %>%
  dplyr::select(-lower, -order_group) %>%
  mutate(across(where(is.numeric), ~ signif(., 5))) %>%
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE))) %>%
  dplyr::select(-c("subsample", "analysis"))

#############################################################
# T) SFS analysis:Compare predicted sel and neu
#############################################################

SFS_upd <- SFS %>%
  dplyr::filter(esm_interval != "[-12.8, 11.3]") %>%
  dplyr::mutate(Spectra = recode_factor(Spectra, observed_neu = "SFS (neutral)", predicted_sel = "SFS (Selected)"))

SFS_proportion_plot <- ggplot(
  SFS_upd %>% dplyr::filter(Spectra != "observed_sel"),
  aes(x = interval, y = Proportion, fill = Spectra)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_bw() +
  facet_wrap(~esm_interval, ncol = 5, nrow = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  labs(x = "Allele frequency", y = "Proportion", fill = "Spectra") +
  ylim(0, 1)

SFS_count_plot <- ggplot(
  SFS_upd %>% dplyr::filter(Spectra != "observed_sel"),
  aes(x = interval, y = Count, fill = Spectra)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_bw() +
  facet_wrap(~esm_interval, ncol = 5, nrow = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  labs(x = "Allele frequency", y = "count", fill = "Spectra")

# Combined
SFS_plots <- ggarrange(SFS_proportion_plot, SFS_count_plot,
  ncol = 1, nrow = 2,
  labels = c("A", "B"), common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(current_path, "spectra/SFS_predicted_combined_", submission, ".png"), plot = SFS_plots, width = 10, height = 8, dpi = 500)

#############################################################

SFS_upd <- SFS %>%
  dplyr::mutate(Spectra = recode_factor(Spectra, observed_neu = "SFS (neutral)", predicted_sel = "SFS (selected)"))

SFS_proportion_plot <- ggplot(
  SFS_upd %>% dplyr::filter(Spectra != "observed_sel"),
  aes(x = interval, y = Proportion, fill = Spectra)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_bw() +
  facet_wrap(~esm_interval, ncol = 6, nrow = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  labs(x = "Allele frequency", y = "Proportion", fill = "Spectra") +
  ylim(0, 1)

SFS_count_plot <- ggplot(
  SFS_upd %>% dplyr::filter(Spectra != "observed_sel"),
  aes(x = interval, y = Count, fill = Spectra)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_bw() +
  facet_wrap(~esm_interval, ncol = 6, nrow = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank()
  ) +
  labs(x = "Allele frequency", y = "count", fill = "Spectra")

SFS_plots <- ggarrange(SFS_proportion_plot, SFS_count_plot,
  ncol = 1, nrow = 2,
  labels = c("A", "B"), common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(current_path, "spectra/SFS_predicted_combined_all_", submission, ".png"), plot = SFS_plots, width = 10, height = 8, dpi = 500)

#############################################################
# T) SFS analysis: Prepare figures for manuscript
#############################################################

fill_colors <- c(
  "[-12.8, -4.8)" = "#440154", "[-4.8, -3.8)" = "#482878", "[-3.8, -3.1)" = "#3E4989", "[-3.1, -2.2)" = "#31688E",
  "[-2.2, -0.9)" = "#26828E", "[-0.9, 0.9)" = "#1F9E89", "[0.9, 2.1)" = "#35B779", "[2.1, 3)" = "#6DCD59",
  "[3, 3.9)" = "#B4DE2C", "[3.9, 11.3]" = "#FDE725", "[3.9, 11.3)" = "#FDE725", "[-12.8, 11.3]" = "grey"
)

DFE_discrete_plot <- ggplot(
  DFE_discrete,
  aes(x = interval, y = proportion, fill = esm_interval)
) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = proportion - LC, ymax = proportion + UC), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = fill_colors) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title.position = "top", nrow = 1)) +
  labs(x = "Scaled selection coefficients", y = "DFE", fill = "Evolutionary score (interval)")

ggsave(paste0(current_path, "DFE_discrete_", submission, "_allsites.png"), plot = DFE_discrete_plot, width = 12, height = 4, dpi = 500)


fill_colors_sub <- fill_colors[!(names(fill_colors) %in% "[-12.8, 11.3]")]

DFE_discrete_upd_plot <- ggplot(
  DFE_discrete %>% dplyr::filter(esm_interval != "[-12.8, 11.3]"),
  aes(x = interval, y = proportion, fill = esm_interval)
) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = proportion - LC, ymax = proportion + UC), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = fill_colors_sub) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title.position = "top", nrow = 2)) +
  labs(x = "Scaled selection coefficients", y = "DFE", fill = "Evolutionary score (interval)        ")
ggsave(paste0(current_path, "DFE_discrete_", submission, "_final.png"), plot = DFE_discrete_upd_plot, width = 12, height = 4, dpi = 500)

#
SFS_DFE <- list()

for (current_range in bin_order) {
  print(current_range)

  current_SFS <- SFS_sel %>% dplyr::filter(esm_interval == current_range)

  SFS_plot <- ggplot(current_SFS, aes(x = interval, y = Proportion, fill = esm_interval)) +
    geom_bar(stat = "identity", alpha = 0.75) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 10)
    ) +
    scale_fill_manual(values = fill_colors) +
    labs(x = "Allele frequency", y = bquote(P[S] * ", " * .(current_range))) +
    ylim(0, 1) +
    theme(legend.position = "none")

  ggsave(paste0(current_path, "marginal/SFS_predicted_", current_range, "_", submission, ".png"), plot = SFS_plot, width = 4, height = 3, dpi = 500)

  current_DFE <- DFE_discrete %>% dplyr::filter(esm_interval == current_range)

  DFE_plot <- ggplot(current_DFE, aes(x = interval, y = proportion, fill = esm_interval)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = proportion - LC, ymax = proportion + UC), width = 0.2, position = position_dodge(0.9)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      plot.background = element_blank()
    ) +
    scale_fill_manual(values = fill_colors) +
    labs(x = "Scaled selection coefficients", y = "DFE") +
    ylim(0, 1) +
    theme(legend.position = "none")

  ggsave(paste0(current_path, "marginal/DFE_predicted_", current_range, "_", submission, ".png"), plot = DFE_plot, width = 4, height = 2.75, dpi = 500)

  DFE_plot <- DFE_plot + labs(x = "", y = "DFE")

  SFS_DFE[[current_range]] <- SFS_plot + inset_element(DFE_plot, left = 0.10, bottom = 0.25, right = 0.80, top = 0.99)

  ggsave(paste0(current_path, "marginal/SFS_DFE_predicted_", current_range, "_", submission, ".png"),
    plot = SFS_DFE[[current_range]], width = 6, height = 3.5, dpi = 500
  )
}

#############################################################

DA_quantiles <- quantile(Scores$esm_score_DA, probs = seq(0, 1, by = 0.1)) %>% unname()
label_DA_df <- data.frame(x = DA_quantiles, y = 700 * 1.05, label = round(DA_quantiles, 1))

Scores_upd <- Scores %>% dplyr::mutate(esm_interval_DA = recode(esm_interval_DA, "[3.9, 11.3)" = "[3.9, 11.3]"))

Histogram_DA_plot <- ggplot(
  Scores_upd %>% dplyr::filter(esm_score_DA != "[-12.8, 11.3]"),
  aes(x = esm_score_DA, fill = esm_interval_DA)
) +
  geom_histogram(alpha = 1, bins = 150) +
  geom_vline(xintercept = DA_quantiles, linetype = "dashed", color = "black") +
  scale_fill_manual(values = fill_colors_sub) +
  geom_text(data = label_DA_df, aes(x = x, y = y, label = label), angle = 90, vjust = -0.5, size = 4, inherit.aes = FALSE) +
  labs(x = "Evolutionary scores of derived alleles", y = "") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0, 800)

# Combined figure
Full_SFS_plot <- ggpubr::ggarrange(SFS_DFE[["[-12.8, 11.3]"]], labels = c("A"), nrow = 1, ncol = 1)

Combined_SFS_plot <- ggpubr::ggarrange(SFS_DFE[["[-12.8, -4.8)"]], SFS_DFE[["[-4.8, -3.8)"]],
  SFS_DFE[["[-3.8, -3.1)"]], SFS_DFE[["[-3.1, -2.2)"]],
  SFS_DFE[["[-2.2, -0.9)"]], SFS_DFE[["[-0.9, 0.9)"]],
  SFS_DFE[["[0.9, 2.1)"]], SFS_DFE[["[2.1, 3)"]],
  SFS_DFE[["[3, 3.9)"]], SFS_DFE[["[3.9, 11.3)"]],
  labels = c("B", "C", "D", "E", "F", "G", "H", "I", "J", "K"), nrow = 5, ncol = 2, common.legend = F
)

Full_DFE_plot <- ggpubr::ggarrange(DFE_discrete_plot, labels = c("L"), nrow = 1, ncol = 1)

LD_decay_plot_upd <- LD_decay_plot + theme(legend.position = "none")

# Main figure
Final_DFE_plot_ext <- ggpubr::ggarrange(Histogram_DA_plot, LD_decay_plot_upd, DFE_discrete_upd_plot,
  nrow = 3, ncol = 1,
  labels = c("A", "B", "C"), heights = c(2, 2, 2.75)
)

ggsave(paste0(path_results, "figures/SFS/DFE_discrete_combined_main_ext_", submission, ".png"),
  Final_DFE_plot_ext,
  width = 10.5, height = 10, dpi = 500
)

# Supplementary figure
Final_SFS_plot <- ggpubr::ggarrange(Full_SFS_plot, Combined_SFS_plot, Full_DFE_plot, nrow = 3, ncol = 1, heights = c(2.5, 10, 2.5))
ggsave(paste0(path_results, "figures/SFS/DFE_discrete_combined_supp_", submission, ".png"), Final_SFS_plot, width = 12, height = 18, dpi = 500)

#############################################################
# SFS analysis; prepare figures for presentation
#############################################################

Scores_final <- Scores_upd %>% dplyr::filter(esm_score_DA != "[-12.8, 11.3]")

for (current_range in unique(Scores_final$esm_interval_DA)) {
  current_col <- fill_colors_sub
  current_col[names(current_col) != current_range] <- "#D3D3D3"

  Histogram_DA_plot <- ggplot(
    Scores_final,
    aes(x = esm_score_DA, fill = esm_interval_DA)
  ) +
    geom_histogram(alpha = 1, bins = 150) +
    geom_vline(xintercept = DA_quantiles, linetype = "dashed", color = "black") +
    scale_fill_manual(values = current_col) +
    geom_text(data = label_DA_df, aes(x = x, y = y, label = label), angle = 90, vjust = -0.5, size = 4, inherit.aes = FALSE) +
    labs(x = "Evolutionary scores (derived alleles)", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    ylim(0, 800)

  ggsave(paste0(current_path, "presentation/Histogram_scores_DA_", current_range, "_", submission, ".png"),
    Histogram_DA_plot,
    width = 12, height = 4, dpi = current_dpi
  )
}

#############################################################
# SFS analysis; Joint SFS
#############################################################

SFS_DFE_joint_sub <- SFS_DFE_joint %>% distinct(shared_param, llik_pval)

#############################################################
# U) Allele frequency plot
#############################################################

input_data <- Missense_list[["esm_scores"]] %>%
  dplyr::mutate(identifier = paste(CHROM_POS, TRANSCRIPT_ID, sep = "_")) %>%
  dplyr::filter(!is.na(esm_score) & !is.na(sift_score) & !is.na(SAP_frequency)) %>%
  dplyr::mutate(
    esm_score_scaled = (esm_score - min(esm_score)) / (max(esm_score) - min(esm_score)),
    sift_score_scaled = (sift_score - min(sift_score)) / (max(sift_score) - min(sift_score))
  ) %>%
  dplyr::select(-c("esm_score", "sift_score")) %>%
  pivot_longer(cols = c(esm_score_scaled, sift_score_scaled), names_to = "evolutionary_score", values_to = "score") %>%
  dplyr::mutate(evolutionary_score = recode(evolutionary_score,
    "esm_score_scaled" = "ESM", "sift_score_scaled" = "SIFT"
  )) %>%
  dplyr::select(c("evolutionary_score", "score", "SAP_frequency", "CHROM"))

output <- data.frame()
coef_df <- data.frame()

for (metric in unique(input_data$evolutionary_score)) {
  input <- input_data %>% dplyr::filter(evolutionary_score == metric, !is.na(score))

  fit <- gam(SAP_frequency ~ CHROM - 1 + s(score, bs = "cr"), data = input)

  temp_coef_df <- summary(fit)

  temp_coef <- data.frame("R2_adj" = temp_coef_df$r.sq, "Deviance" = temp_coef_df$dev.expl, "chi_sq" = temp_coef_df$chi.sq) %>%
    dplyr::mutate("evolutionary_score" = metric) %>%
    dplyr::mutate_if(is.numeric, round, digits = 3)

  coef_df <- rbind(temp_coef, coef_df)

  fit_coef <- summary(fit)$p.coeff
  fit_se <- summary(fit)$se

  temp_output <- predict_gam(fit, values = list(x = seq(min(input$score, na.rm = TRUE), max(input$score, na.rm = TRUE), 0.001), CHROM = "1")) %>%
    dplyr::select(score, fit, se.fit) %>%
    as.data.frame() %>%
    dplyr::mutate("evolutionary_score" = metric)

  output <- rbind(output, temp_output)
}

#############################################################

R2_sift <- coef_df %>%
  dplyr::filter(evolutionary_score == "SIFT") %>%
  pull(R2_adj)
R2_esm <- coef_df %>%
  dplyr::filter(evolutionary_score == "ESM") %>%
  pull(R2_adj)

GAM_v1 <- ggplot(
  output,
  aes(x = score, y = fit, color = evolutionary_score, fill = evolutionary_score)
) +
  geom_smooth_ci(evolutionary_score, size = 1, linetype = "solid", ci_z = qnorm(0.975), ci_alpha = 0.2) +
  theme_bw() +
  labs(x = "Evolutionary scores of ALT alleles", y = "Allele frequency in diversity panel", color = "score", fill = "score") +
  theme(
    legend.position = "left",
    legend.key.size = unit(0.75, "cm"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white")
  ) +
  geom_text(aes(x = 0.4, y = 0.95, label = paste("ESM: R2.adj=", R2_esm, "and SIFT: R2.adj =", R2_sift)),
    inherit.aes = FALSE, hjust = 0, size = 3
  ) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  scale_color_viridis_d(option = "D") +
  ylim(0, 1)

allele_freq_hist_v1 <- ggplot(
  input_data %>% dplyr::filter(evolutionary_score == "ESM"),
  aes(x = SAP_frequency, fill = evolutionary_score, color = evolutionary_score)
) +
  geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
  geom_density(alpha = 0.25, color = NA) +
  scale_fill_manual(values = "grey", name = "Score") +
  scale_color_manual(values = "grey", name = "Score") +
  theme_bw() +
  labs(x = "", y = "", fill = "Score") +
  coord_flip() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_flip(expand = FALSE) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "Distribution allele frequency")

esm_freq_hist_v1 <- ggplot(
  input_data,
  aes(x = score, fill = evolutionary_score, color = evolutionary_score)
) +
  geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
  facet_grid(~evolutionary_score) +
  scale_fill_viridis_d(option = "D", name = "Score") +
  scale_color_viridis_d(option = "D", name = "Score") +
  theme_bw() +
  labs(x = NULL, y = "Distribution of evolutionary scores") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )

main_v1 <- GAM_v1 + theme(plot.margin = margin(0, 120, 0, 0))

combined_v1 <- esm_freq_hist_v1 / main_v1 + plot_layout(heights = c(1.5, 5))

combined_v1_plot <- combined_v1 + inset_element(allele_freq_hist_v1, left = 1.01, right = 1.30, bottom = -0.01, top = 1, clip = FALSE)
ggsave(paste0(path_results, "figures/allele_frequency/AF_combined_ALT_", submission, ".png"),
  combined_v1_plot,
  width = 7, height = 5.5, dpi = 500
)

#############################################################
# R) Variance Inflation analysis
#############################################################

cols <- c(
  "[-12.8, -4.8)" = "#440154", "[-4.8, -3.8)" = "#482878", "[-3.8, -3.1)" = "#3E4989", "[-3.1, -2.2)" = "#31688E", "[-2.2, -0.9)" = "#26828E",
  "[-0.9, 0.9)" = "#1F9E89", "[0.9, 2.1)" = "#35B779", "[2.1, 3)" = "#6DCD59", "[3, 3.9)" = "#B4DE2C", "[3.9, 11.3)" = "#FDE725",
  "[3.9, 11.3]" = "#FDE725", "PC"
)

loads <- readRDS(paste0(path_input, "kinships/fixed_effects/GW_load_DA.rds"))

phenotypes <- readRDS(paste0(path_input, "prepared_data/Phenotypes.rds"))

info <- loads %>% distinct(esm_interval, q_interval)

current_path <- paste0(paste0(path, "results/"), "figures/final_figures/")

VIF_PCA_df <- data.frame()

for (current_interval in unique(loads$q_interval)) {
  loads_sub <- loads %>%
    dplyr::filter(q_interval == current_interval) %>%
    left_join(phenotypes[[1]], by = "id")

  dummy_y <- rnorm(nrow(loads_sub))

  model <- lm(dummy_y ~ G + P + PC1 + PC2 + PC3, data = loads_sub)

  temp <- car::vif(model) %>%
    as.data.frame() %>%
    dplyr::mutate("q_interval" = current_interval) %>%
    rownames_to_column("covariate") %>%
    rename("VIF" = ".")

  VIF_PCA_df <- rbind(temp, VIF_PCA_df)
}

VIF_PCA_df_upd <- VIF_PCA_df %>%
  left_join(info, by = "q_interval") %>%
  mutate(esm_interval_upd = ifelse(covariate %in% c("PC1", "PC2", "PC3"), "PC", esm_interval)) %>%
  mutate(esm_interval = factor(esm_interval, levels = names(cols)))

S4_plot <- ggplot(VIF_PCA_df_upd, aes(x = covariate, y = VIF, fill = esm_interval_upd)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  labs(x = "Covariates", y = "Variance inflation factor") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  facet_wrap(~esm_interval, nrow = 2, ncol = 5, scales = "free_x")

###########################################################
### Final results
#############################################################

current_path <- paste0(path_results, "figures/final_figures/")

###########################################################
# Final figures
#############################################################

# Main figures
# M1)
M1_plot <- ggarrange(PCA_cluster_plot_upd, combined_v1_plot, ncol = 2, labels = c("A", "B"), widths = c(1, 1.25))
ggsave(paste0(current_path, "M12_", submission, ".png"), M1_plot, width = 14, height = 5, dpi = 500)

# M2)
M2_plot <- Final_DFE_plot_ext
ggsave(paste0(current_path, "M3_", submission, ".png"), M2_plot, width = 8, height = 10, dpi = 500)

# M3)-see other script

# M4)-see other script

### Supplementary figures

# S1)
S1_plot <- S1_combined
ggsave(paste0(current_path, "S1_", submission, ".png"), S1_plot, width = 12, height = 9, dpi = 500)

# S2)
S2_plot <- Final_SFS_plot
ggsave(paste0(current_path, "S2_", submission, ".png"), S2_plot, width = 12, height = 18, dpi = 500)

# S3)
S3_plot <- plots
ggsave(paste0(current_path, "S3_", submission, ".png"), S3_plot, height = 20, width = 10, dpi = 500)

# S4)
S4_plot <- S4_plot
ggsave(paste0(current_path, "S4_", submission, ".png"), S4_plot, width = 10, height = 5, dpi = 500)

###########################################################
# Final tables
#############################################################

# T0)
T0_table <- LD_anova
ggsave(paste0(current_path, "T0.pdf"), plot = T0_table, width = 8, height = 3)

# T1)
T1_table <- SFS_param %>% tableGrob(rows = NULL)
ggsave(paste0(current_path, "T1.pdf"), plot = T1_table, width = 12, height = 4)

# T2)
T2_table <- Mean_part_table
ggsave(paste0(current_path, "T2.pdf"), plot = T2_table, width = 9, height = 3)

# T3)
T3_table <- Var_part_table
ggsave(paste0(current_path, "T3.pdf"), plot = T3_table, width = 9, height = 3)

# T4)
T4_table <- M0_PA %>% tableGrob(rows = NULL)
ggsave(paste0(current_path, "T4.pdf"), plot = T4_table, width = 6, height = 4)

# T5)
T5_table <- Val_part_table
ggsave(paste0(current_path, "T5.pdf"), plot = T5_table, width = 13, height = 4)

#############################################################
#############################################################
