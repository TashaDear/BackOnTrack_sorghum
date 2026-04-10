#############################################################
library("data.table")
library("dplyr")
library("SNPRelate")
library("foreach")
library("doParallel")
library("foreach")
library("doParallel")
registerDoParallel(30)
#############################################################
# A) Paths
#############################################################

set.seed(0)

submission <- "submission1"

path <- "/home/tasha/backontrack/"
path_input <- paste0(path, "input/")
path_output <- paste0(path, "output/")
path_results <- paste0(path, "results/")

#############################################################
# B) Functions
#############################################################

invSqrt <- function(mat) {
  stopifnot(isSymmetric(mat))

  ED <- eigen(mat)
  L <- ED$vectors %*% diag(1 / sqrt(ED$values)) %*% t(ED$vectors)

  return(L)
}

#############################################################
# C) Prepare matrices
#############################################################

Geno <- readRDS(paste0(path_input, "kinships/Geno_list.rds"))[["Geno"]]

Q <- prcomp(Geno, center = T, scale = F)$x[, 1:3]
Q <- cbind(Intercept = 1, Q)
G <- readRDS(paste0(path_input, "kinships/random_effects/GRMs.rds"))[["Geno"]]
genotypes <- rownames(Q)
n <- length(genotypes)
G <- G[genotypes, genotypes]
diag(G) <- diag(G) + 1e-5

# Projection matrices
P <- list()
P[["Int"]] <- diag(n) - matrix(1, nrow = n, ncol = n) / n
P[["Q"]] <- diag(n) - Q %*% solve(crossprod(Q)) %*% t(Q)

Ginv <- solve(G)
H <- Q %*% solve(t(Q) %*% Ginv %*% Q) %*% t(Q) %*% Ginv
P[["Q+K"]] <- invSqrt(G) %*% (diag(n) - H)

max_dist <- 100000

#############################################################
# D) Load input
#############################################################

pair_1 <- readRDS(paste0(path_input, "/Missense_list.rds"))[["Missense_AA_df"]] %>%
  dplyr::filter(Degeneracy == "0" & !is.na(esm_score))

pair_2 <- readRDS(paste0(path_input, "LDdecay/Input_LD_list.rds"))[["SNPs_within"]] %>%
  rename("CHROM" = "chromosome") %>%
  sample_frac(0.25)

#############################################################
# E) LD across samples
#############################################################

output <- list()

for (chrom in unique(pair_1$CHROM)) {
  print(chrom)

  pair_1_upd <- pair_1 %>% filter(CHROM == chrom)
  pair_2_upd <- pair_2 %>% filter(CHROM == chrom)

  X_raw <- Geno[, unique(c(pair_1_upd$CHROM_POS, pair_2_upd$CHROM_POS))]

  temp_results <- foreach(i = pair_1_upd$CHROM_POS, .combine = "rbind", .packages = "SNPRelate") %dopar% {
    pos_i <- pair_1_upd %>%
      filter(CHROM_POS == i) %>%
      pull(POS) %>%
      as.numeric()

    nearby_snps <- pair_2_upd %>% filter(position >= pos_i - max_dist, position <= pos_i + max_dist, CHROM_POS != i)

    ld_list <- lapply(names(P), function(adj) {
      X_adj <- P[[adj]] %*% X_raw
      X_adj <- scale(X_adj, center = FALSE, scale = sqrt(apply(X_adj, 2, crossprod)))

      temp <- lapply(nearby_snps$CHROM_POS, function(j) {
        pos_j <- nearby_snps %>%
          filter(CHROM_POS == j) %>%
          pull(position)
        r <- sum(X_adj[, i] * X_adj[, j])
        r <- min(max(r, -1), 1)

        data.frame(adjustment = adj, CHROM = chrom, CHROM_POS = i, snp = i, D = abs(pos_i - pos_j), r = r, r2 = r^2)
      })

      bind_rows(temp)
    })

    bind_rows(ld_list)
  }

  temp_split <- split(temp_results, temp_results$adjustment)

  for (adj in names(temp_split)) {
    output[[adj]] <- bind_rows(output[[adj]], temp_split[[adj]])
  }
}

saveRDS(output, paste0(path_input, "LDdecay/LD_prioritized_", submission, ".rds"))

#############################################################
