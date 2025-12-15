# R/modules_stats.R
# ==============================================================================
# STATISTICAL INFERENCE ENGINE
# Description: Pure functions for Robust Partial Correlation, Bootstrap & Permutation.
# Dependencies: corpcor, dplyr, foreach
# ==============================================================================

library(corpcor)
library(dplyr)
library(foreach)

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