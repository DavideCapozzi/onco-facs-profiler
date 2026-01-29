# R/modules_network.R
# ==============================================================================
# NETWORK INFERENCE MODULE
# Description: Partial Correlation (Shrinkage), Bootstrap, and Topology.
# Dependencies: corpcor, foreach, igraph, dplyr
# ==============================================================================

library(corpcor)
library(dplyr)
library(foreach)

#' @title Infer Robust Partial Correlation Network
#' @description 
#' Uses Schafer-Strimmer Shrinkage to directly estimate partial correlations.
#' 
#' @param mat Numeric matrix (Samples x Features).
#' @param fixed_lambda Optional numeric. If provided, enforces this shrinkage intensity.
#'        If NULL, lambda is estimated analytically (James-Stein).
#' @return Partial correlation matrix with diagonal = 1.
infer_network_pcor <- function(mat, fixed_lambda = NULL) {
  
  # Direct calculation using corpcor's optimized function
  # This performs James-Stein shrinkage on covariance and efficient inversion
  if (is.null(fixed_lambda)) {
    # Default: Estimate lambda from data
    pcor_mat <- corpcor::pcor.shrink(mat, verbose = FALSE)
  } else {
    # Enforce fixed lambda (e.g., for consistent bootstrap)
    pcor_mat <- corpcor::pcor.shrink(mat, lambda = fixed_lambda, verbose = FALSE)
  }
  
  return(as.matrix(pcor_mat))
}

#' @title Single Bootstrap Worker
#' @description 
#' Worker function for parallel loops. Supports fixed lambda injection.
#' 
#' @param mat Numeric matrix (Samples x Features).
#' @param n Number of samples to resample.
#' @param lambda_val Optional fixed lambda value to pass to inference.
boot_worker_pcor <- function(mat, n, lambda_val = NULL) {
  idx <- sample(1:n, n, replace = TRUE)
  boot_sample <- mat[idx, ]
  
  # Safety Check: Zero variance (constant columns) causes correlation failure
  # If any column is constant in this resample, discard the iteration.
  col_vars <- apply(boot_sample, 2, var, na.rm = TRUE)
  if (any(col_vars == 0 | is.na(col_vars))) {
    return(NULL)
  }
  
  tryCatch({
    return(infer_network_pcor(boot_sample, fixed_lambda = lambda_val))
  }, error = function(e) return(NULL))
}

#' @title Check Bootstrap Yield
#' @description Warnings if too many bootstrap iterations were discarded due to zero variance.
#' @param res_obj The result list from aggregate_boot_results.
#' @param n_req The original number of requested bootstraps.
#' @param label String label for the group (e.g., "Control").
check_boot_yield <- function(res_obj, n_req, label) {
  if (is.null(res_obj) || is.null(res_obj$n_boot_valid)) return()
  
  if (res_obj$n_boot_valid < (n_req * 0.8)) {
    warning(sprintf(
      "[WARN] Low bootstrap yield for %s: %d/%d valid iterations. Results may be unstable.", 
      label, res_obj$n_boot_valid, n_req
    ))
  } else {
    message(sprintf("   -> %s Bootstrap Yield: %d/%d (OK)", label, res_obj$n_boot_valid, n_req))
  }
}

#' @title Aggregate Bootstrap Results
#' @description Calculates edge stability and mean weights based on CI.
aggregate_boot_results <- function(boot_list, alpha = 0.05) {
  valid_boots <- Filter(Negate(is.null), boot_list)
  n_boots <- length(valid_boots)
  if (n_boots == 0) stop("All bootstrap iterations failed.")
  
  p <- nrow(valid_boots[[1]])
  nodes <- rownames(valid_boots[[1]])
  
  adj_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  mean_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  stab_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      vals <- sapply(valid_boots, function(m) m[i, j])
      ci <- quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      
      # Select edge if CI does not cross zero
      if ((ci[1] > 0) || (ci[2] < 0)) {
        adj_mat[i, j] <- 1; adj_mat[j, i] <- 1
        mu <- mean(vals, na.rm = TRUE)
        mean_mat[i, j] <- mu; mean_mat[j, i] <- mu
      }
      stab_mat[i, j] <- sum(sign(vals) == sign(mean(vals))) / n_boots
      stab_mat[j, i] <- stab_mat[i, j]
    }
  }
  return(list(adj = adj_mat, weights = mean_mat, stability = stab_mat, n_boot_valid = n_boots))
}

#' @title Calculate Topology Metrics
get_topology_metrics <- function(adj_mat, weight_mat) {
  requireNamespace("igraph", quietly = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  if (igraph::vcount(g) == 0) return(NULL)
  
  data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )
}

#' @title Export Edges for Cytoscape
export_cytoscape_edges <- function(adj_mat, weight_mat) {
  requireNamespace("igraph", quietly = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  edge_list <- igraph::as_data_frame(g, what = "edges")
  
  if (nrow(edge_list) == 0) return(data.frame())
  
  edge_list %>%
    rowwise() %>%
    mutate(
      Weight = weight_mat[from, to],
      Interaction = ifelse(Weight > 0, "Co-occurrence", "Mutual-Exclusion"),
      Abs_Weight = abs(Weight)
    ) %>%
    dplyr::rename(Source = from, Target = to) %>%
    ungroup()
}

#' @title Infer Differential Network
#' @description 
#' Calculates the difference in correlation structure between two conditions.
#' Returns a table of edges that change significantly between Case and Control.
#' 
#' @param mat_ctrl Numeric matrix for Control group.
#' @param mat_case Numeric matrix for Case group.
#' @param method Correlation method ("spearman" or "pearson").
#' @param n_bootstrap Number of bootstraps (placeholder for future implementation).
#' @return A list containing the edge table.
infer_differential_network <- function(mat_ctrl, mat_case, method = "spearman", n_bootstrap = 100) {
  
  # Safety Checks
  if (ncol(mat_ctrl) != ncol(mat_case)) stop("Matrices must have same columns.")
  
  # Compute Correlations
  # Use use="pairwise.complete.obs" to handle potential NAs robustly
  cor_ctrl <- cor(mat_ctrl, method = method, use = "pairwise.complete.obs")
  cor_case <- cor(mat_case, method = method, use = "pairwise.complete.obs")
  
  # Replace NAs with 0 (if variance was 0 for some features)
  cor_ctrl[is.na(cor_ctrl)] <- 0
  cor_case[is.na(cor_case)] <- 0
  
  # Calculate Difference matrix
  diff_mat <- cor_case - cor_ctrl
  
  # Flatten to Edge List (Upper triangle only)
  edges <- list()
  markers <- colnames(mat_ctrl)
  
  count <- 1
  for (i in 1:(length(markers)-1)) {
    for (j in (i+1):length(markers)) {
      m1 <- markers[i]
      m2 <- markers[j]
      
      w_ctrl <- cor_ctrl[i, j]
      w_case <- cor_case[i, j]
      delta  <- diff_mat[i, j]
      
      # Filter noise: Keep edge only if correlation is strong in at least one condition
      # OR if the shift is massive (> 0.4)
      if (abs(w_ctrl) > 0.3 || abs(w_case) > 0.3 || abs(delta) > 0.4) {
        edges[[count]] <- data.frame(
          Node1 = m1,
          Node2 = m2,
          Weight_Ctrl = round(w_ctrl, 3),
          Weight_Case = round(w_case, 3),
          Diff_Score = round(delta, 3),
          Type = ifelse(delta > 0, "Strengthened", "Weakened"),
          stringsAsFactors = FALSE
        )
        count <- count + 1
      }
    }
  }
  
  if (length(edges) > 0) {
    edges_df <- do.call(rbind, edges)
    # Sort by absolute difference
    edges_df <- edges_df[order(-abs(edges_df$Diff_Score)), ]
  } else {
    edges_df <- NULL
  }
  
  return(list(diff_edges = edges_df))
}