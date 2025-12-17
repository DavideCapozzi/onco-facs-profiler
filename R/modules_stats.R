# R/modules_stats.R
# ==============================================================================
# STATISTICAL INFERENCE ENGINE
# Description: Pure functions for Robust Partial Correlation, Bootstrap & Permutation.
# Dependencies: corpcor, dplyr, foreach
# ==============================================================================

library(corpcor)
library(dplyr)
library(foreach)

#' @title Assess MANOVA Assumptions
#' @description 
#' Checks Multivariate Normality (via Mahalanobis QQ) and 
#' Homogeneity of Multivariate Dispersion (via Betadisper).
#' 
#' @param ilr_data Dataframe with metadata and ILR coordinates.
#' @param group_col String name of the grouping column.
#' @param metadata_cols Vector of metadata column names to exclude.
#' @return A list containing statistical test results and metrics.
check_manova_assumptions <- function(ilr_data, group_col = "Group", metadata_cols = c("Patient_ID", "Group")) {
  
  requireNamespace("vegan", quietly = TRUE)
  
  # 1. Prepare Data
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups <- as.factor(ilr_data[[group_col]])
  
  results <- list()
  
  # 2. Multivariate Normality (Proxy via Mahalanobis Distance)
  # Theoretical quantiles vs Squared Mahalanobis distances
  # If data is MVN, squared Mahalanobis distances follow Chi-Square(df = n_features)
  
  mu <- colMeans(ilr_mat)
  cov_mat <- cov(ilr_mat)
  
  # Handle singular covariance (n < p) using generalized inverse if needed
  # But for Mahalanobis, we usually need n > p. If n < p, we skip strict MVN test.
  if(nrow(ilr_mat) > ncol(ilr_mat)) {
    d2 <- mahalanobis(ilr_mat, mu, cov_mat)
    
    # Shapiro test on the Mahalanobis distances (Are they Chi-Sq distributed?)
    # Note: This is an approximation. 
    # Technically, we test if d2 follows Chi-sq. A simple proxy is correlation.
    results$mahalanobis_dist <- d2
    results$mvn_proxy_cor <- cor(sort(d2), qchisq(ppoints(nrow(ilr_mat)), df = ncol(ilr_mat)))
  } else {
    results$mahalanobis_dist <- NULL
    results$mvn_proxy_cor <- NA
    warning("[Stats] N < P detected. Skipping strict Mahalanobis-based MVN check.")
  }
  
  # 3. Homogeneity of Multivariate Dispersions (PERMDISP)
  # This is the robust alternative to Box's M
  # Measures distance of each point to its group centroid
  dist_mat <- dist(ilr_mat, method = "euclidean")
  mod_disp <- vegan::betadisper(dist_mat, groups)
  
  # ANOVA on the dispersions (Do groups have different "tightness"?)
  test_disp <- anova(mod_disp)
  
  results$homogeneity_pval <- test_disp$`Pr(>F)`[1]
  results$dispersions <- mod_disp$distances
  results$groups <- groups
  
  return(results)
}

#' @title Run PERMANOVA (Robust Alternative)
#' @description 
#' Performs Permutational Multivariate Analysis of Variance (adonis2).
#' Valid when Normality assumptions are violated.
#' 
#' @param ilr_data Dataframe with metadata and ILR coordinates.
#' @param group_col Grouping column.
#' @param n_perm Number of permutations.
#' @return A list with the adonis table and p-value.
run_coda_permanova <- function(ilr_data, group_col = "Group", n_perm = 999) {
  requireNamespace("vegan", quietly = TRUE)
  
  # Prepare formula
  # We construct the formula dynamically: ilr_mat ~ Group
  metadata_cols <- c("Patient_ID", group_col)
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups <- ilr_data[[group_col]]
  
  message(sprintf("   [Stats] Running PERMANOVA (adonis2) with %d permutations...", n_perm))
  
  # Run adonis2 (Euclidean on ILR = Aitchison on Simplex)
  # by="margin" tests marginal effects (Type III SS equivalent)
  res_adonis <- vegan::adonis2(ilr_mat ~ groups, method = "euclidean", permutations = n_perm)
  
  return(res_adonis)
}

#' @title Run CoDa MANOVA (Multivariate Analysis of Variance)
#' @description 
#' Performs MANOVA on ILR-transformed data to test global group differences.
#' Calculates Canonical Variate for visualization and correlations with CLR 
#' markers for interpretability (Loadings).
#' 
#' @param ilr_data Dataframe containing metadata ('Group') and ILR coordinates.
#' @param clr_data Dataframe containing metadata and CLR markers (for interpretation).
#' @param metadata_cols Vector of column names to exclude from the matrix.
#' @return A list containing the MANOVA object, summary, viz scores, and marker loadings.
run_coda_manova <- function(ilr_data, clr_data, metadata_cols = c("Patient_ID", "Group")) {
  
  # 1. Prepare ILR Data (For Statistics)
  ilr_mat <- as.matrix(ilr_data[, !names(ilr_data) %in% metadata_cols])
  groups  <- as.factor(ilr_data$Group)
  
  # 2. Run MANOVA (on ILR)
  manova_res <- manova(ilr_mat ~ groups)
  summ_res <- summary(manova_res, test = "Pillai")
  
  stats_df <- as.data.frame(summ_res$stats)
  pval <- stats_df[1, "Pr(>F)"]
  pillai <- stats_df[1, "Pillai"]
  
  # 3. Canonical Discriminant Analysis (LDA)
  # This finds the axis of separation
  requireNamespace("MASS", quietly = TRUE)
  lda_res <- MASS::lda(groups ~ ilr_mat)
  
  # Get scores (The position of each patient on the separation axis)
  scores <- predict(lda_res, as.data.frame(ilr_mat))$x
  
  # 4. Calculate Loadings (Structure Coefficients) - THE EXPLAINABILITY PART
  # We correlate the CLR markers (interpretable) with the LDA Scores (separation axis)
  # Ensure clr_data is aligned and numeric
  clr_mat <- as.matrix(clr_data[, !names(clr_data) %in% metadata_cols])
  
  # Check alignment
  if(nrow(clr_mat) != nrow(scores)) stop("Row mismatch between CLR data and LDA scores.")
  
  # Calculate correlation (Robust loadings)
  # Result is a vector: correlation of each marker with the separation axis
  loadings <- cor(clr_mat, scores[,1], use = "pairwise.complete.obs")
  
  loadings_df <- data.frame(
    Marker = rownames(loadings),
    Loading = as.numeric(loadings)
  ) %>%
    arrange(desc(abs(Loading))) # Sort by importance
  
  # 5. Output
  plot_data <- data.frame(
    Patient_ID = ilr_data$Patient_ID,
    Group = groups,
    Canonical_Variate_1 = scores[,1]
  )
  
  return(list(
    model = manova_res,
    summary = stats_df,
    p_value = pval,
    pillai_stat = pillai,
    plot_data = plot_data,
    loadings = loadings_df # This is the new part
  ))
}

#' @title Convert Precision Matrix to Partial Correlation
#' @description 
#' standardizes the precision matrix (inverse covariance) to partial correlations.
#' Formula: rho_ij = -omega_ij / sqrt(omega_ii * omega_jj)
#' 
#' @param prec_mat A precision matrix (inverse covariance).
#' @return A partial correlation matrix with diagonal = 1.
prec2part <- function(prec_mat) {
  d <- diag(prec_mat)
  # Outer product to get denominator matrix
  denom <- sqrt(outer(d, d))
  # Apply formula
  part <- -prec_mat / denom
  diag(part) <- 1
  return(part)
}

#' @title Robust Estimation of Partial Correlation
#' @description 
#' Uses Schafer-Strimmer Shrinkage to estimate the covariance matrix, 
#' then inverts it to get partial correlations. This is robust for high-dimensional data.
#' 
#' @param mat A numeric matrix (n_samples x n_markers).
#' @return A symmetric partial correlation matrix (p x p).
estimate_pcor_robust <- function(mat) {
  # Estimate covariance with shrinkage (always positive definite)
  # verbose = FALSE to keep logs clean
  cov_shrink <- corpcor::cov.shrink(mat, verbose = FALSE)
  
  # Calculate precision matrix (pseudo-inverse if necessary, though shrink usually invertible)
  prec_mat <- corpcor::invcov.shrink(mat, verbose = FALSE)
  
  # Convert to partial correlation
  pcor_mat <- prec2part(prec_mat)
  
  return(as.matrix(pcor_mat))
}

#' @title Single Bootstrap Iteration (Worker Function)
#' @description 
#' Performs one iteration of resampling and estimation. 
#' Intended to be called inside a parallel loop.
#' 
#' @param mat Input data matrix.
#' @param n Number of samples to draw (usually nrow(mat)).
#' @return A partial correlation matrix for the bootstrap sample.
run_single_bootstrap <- function(mat, n) {
  # Resample with replacement
  idx <- sample(1:n, n, replace = TRUE)
  boot_sample <- mat[idx, ]
  
  # Estimate
  tryCatch({
    return(estimate_pcor_robust(boot_sample))
  }, error = function(e) {
    # Fail gracefully in parallel workers
    return(NULL) 
  })
}

#' @title Aggregate Bootstrap Results
#' @description 
#' Calculates edge stability (probability of non-zero) and mean weights from bootstrap list.
#' 
#' @param boot_list A list of matrices returned by the bootstrap loop.
#' @param alpha Significance level for confidence intervals (default 0.05).
#' @return A list containing:
#'   - adj_matrix: Binary adjacency matrix (stable edges).
#'   - weight_matrix: Mean partial correlation of stable edges.
#'   - stability_matrix: Frequency of edge detection (0-1).
aggregate_boot_results <- function(boot_list, alpha = 0.05) {
  
  # Remove failed iterations (NULLs)
  valid_boots <- Filter(Negate(is.null), boot_list)
  n_boots <- length(valid_boots)
  
  if (n_boots == 0) stop("All bootstrap iterations failed.")
  
  # Get dimensions
  p <- nrow(valid_boots[[1]])
  nodes <- rownames(valid_boots[[1]])
  
  # Initialize output matrices
  adj_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  mean_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  stab_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  
  # Iterate over upper triangle (symmetric)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      
      # Extract vector of correlations for edge i-j across all bootstraps
      vals <- sapply(valid_boots, function(m) m[i, j])
      
      # 1. Stability: How often is the sign consistent?
      # Alternatively, use Quantile intervals (Stability Selection via CI)
      ci <- quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      
      # Check if CI excludes zero (Significance)
      is_significant <- (ci[1] > 0) || (ci[2] < 0)
      
      if (is_significant) {
        adj_mat[i, j] <- 1; adj_mat[j, i] <- 1
        mu <- mean(vals, na.rm = TRUE)
        mean_mat[i, j] <- mu; mean_mat[j, i] <- mu
      }
      
      # Simple stability score (consistency of sign matching the mean)
      # This is just for reporting, the CI determines the edge
      stab_mat[i, j] <- sum(sign(vals) == sign(mean(vals))) / n_boots
      stab_mat[j, i] <- stab_mat[i, j]
    }
  }
  
  return(list(
    adj = adj_mat,
    weights = mean_mat,
    stability = stab_mat,
    n_boot_valid = n_boots
  ))
}

#' @title Calculate Network Topology Metrics
#' @description 
#' Converts an adjacency/weight matrix pair into an igraph object and calculates
#' node-level centrality metrics (Degree, Betweenness, Closeness).
#' 
#' @param adj_mat Binary adjacency matrix.
#' @param weight_mat Partial correlation matrix (symmetric).
#' @return A dataframe of node metrics.
get_topology_metrics <- function(adj_mat, weight_mat) {
  requireNamespace("igraph", quietly = TRUE)
  
  # create graph from adjacency
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # assign weights (absolute value for strength, but centrality usually uses connectivity)
  # We use the adjacency for structure, weights can be added as attributes
  
  if (igraph::vcount(g) == 0) return(NULL)
  
  metrics <- data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )
  
  return(metrics)
}

#' @title Format Network for Cytoscape
#' @description 
#' Extracts edges from the matrices and formats them for Cytoscape import.
#' 
#' @param adj_mat Binary adjacency matrix.
#' @param weight_mat Partial correlation matrix.
#' @return A dataframe with Source, Target, Weight, and Interaction columns.
format_cytoscape_edges <- function(adj_mat, weight_mat) {
  requireNamespace("igraph", quietly = TRUE)
  
  # Use igraph to extract edge list easily
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = NULL, diag = FALSE)
  edge_list <- igraph::as_data_frame(g, what = "edges")
  
  if (nrow(edge_list) == 0) return(data.frame())
  
  # Populate weights and interaction types
  out_df <- edge_list %>%
    rowwise() %>%
    mutate(
      Weight = weight_mat[from, to],
      Interaction = ifelse(Weight > 0, "Co-occurrence", "Mutual-Exclusion"),
      Abs_Weight = abs(Weight)
    ) %>%
    dplyr::rename(Source = from, Target = to) %>%
    ungroup()
  
  return(out_df)
}