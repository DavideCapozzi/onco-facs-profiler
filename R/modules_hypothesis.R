# R/modules_hypothesis.R
# ==============================================================================
# HYPOTHESIS TESTING MODULE
# Description: MANOVA, PERMANOVA, and Assumption Checks for CoDa.
# Dependencies: vegan, MASS, dplyr
# ==============================================================================

library(dplyr)

#' @title Assess MANOVA Assumptions
#' @description Checks Multivariate Normality and Homogeneity of Dispersion.
assess_manova_assumptions <- function(ilr_data, group_col = "Group", metadata_cols = c("Patient_ID", "Group")) {
  requireNamespace("vegan", quietly = TRUE)
  
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups <- as.factor(ilr_data[[group_col]])
  
  results <- list()
  
  # 1. MVN Proxy (Mahalanobis)
  mu <- colMeans(ilr_mat)
  cov_mat <- cov(ilr_mat)
  
  if(nrow(ilr_mat) > ncol(ilr_mat)) {
    d2 <- mahalanobis(ilr_mat, mu, cov_mat)
    results$mahalanobis_dist <- d2
    results$mvn_proxy_cor <- cor(sort(d2), qchisq(ppoints(nrow(ilr_mat)), df = ncol(ilr_mat)))
  } else {
    results$mahalanobis_dist <- NULL
    results$mvn_proxy_cor <- NA
    warning("[Stats] N < P detected. Skipping strict Mahalanobis-based MVN check.")
  }
  
  # 2. Homogeneity (PERMDISP)
  dist_mat <- dist(ilr_mat, method = "euclidean")
  mod_disp <- vegan::betadisper(dist_mat, groups)
  test_disp <- anova(mod_disp)
  
  results$homogeneity_pval <- test_disp$`Pr(>F)`[1]
  results$dispersions <- mod_disp$distances
  results$groups <- groups
  
  return(results)
}

#' @title Run PERMANOVA (Robust)
#' @description Performs Permutational Multivariate Analysis of Variance (adonis2).
test_coda_permanova <- function(ilr_data, group_col = "Group", n_perm = 999) {
  requireNamespace("vegan", quietly = TRUE)
  
  metadata_cols <- c("Patient_ID", group_col)
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups <- ilr_data[[group_col]]
  
  message(sprintf("   [Stats] Running PERMANOVA (adonis2) with %d permutations...", n_perm))
  
  res_adonis <- vegan::adonis2(ilr_mat ~ groups, method = "euclidean", permutations = n_perm)
  return(res_adonis)
}

#' @title Run CoDa MANOVA (Parametric)
#' @description Performs standard MANOVA and calculates LDA loadings for interpretability.
test_coda_manova <- function(ilr_data, clr_data, metadata_cols = c("Patient_ID", "Group")) {
  requireNamespace("MASS", quietly = TRUE)
  
  # Prepare Data
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups  <- as.factor(ilr_data$Group)
  
  # MANOVA
  manova_res <- manova(ilr_mat ~ groups)
  summ_res <- summary(manova_res, test = "Pillai")
  stats_df <- as.data.frame(summ_res$stats)
  
  # LDA for Separation & Loadings
  lda_res <- MASS::lda(groups ~ ilr_mat)
  scores <- predict(lda_res, as.data.frame(ilr_mat))$x
  
  # Loadings (Correlation of CLR markers with LDA axis)
  clr_mat <- as.matrix(clr_data[, !names(clr_data) %in% metadata_cols])
  if(nrow(clr_mat) != nrow(scores)) stop("Row mismatch between CLR data and LDA scores.")
  
  loadings <- cor(clr_mat, scores[,1], use = "pairwise.complete.obs")
  loadings_df <- data.frame(Marker = rownames(loadings), Loading = as.numeric(loadings)) %>%
    arrange(desc(abs(Loading)))
  
  plot_data <- data.frame(Patient_ID = ilr_data$Patient_ID, Group = groups, Canonical_Variate_1 = scores[,1])
  
  return(list(
    model = manova_res,
    summary = stats_df,
    p_value = stats_df[1, "Pr(>F)"],
    pillai_stat = stats_df[1, "Pillai"],
    plot_data = plot_data,
    loadings = loadings_df
  ))
}