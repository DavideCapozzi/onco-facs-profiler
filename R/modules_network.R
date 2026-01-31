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

#' @title Run Robust Differential Network Analysis (Bootstrap + Permutation)
#' @description 
#' Combines Bootstrap (for edge stability) and Permutation (for differential significance).
#' This reproduces the rigorous logic: Filter unstable edges -> Test differences on stable ones.
#' 
#' @param mat_ctrl Reference matrix.
#' @param mat_case Test matrix.
#' @param n_boot Number of bootstrap iterations for stability (default 100).
#' @param n_perm Number of permutations for significance (default 1000).
#' @param seed Random seed.
#' @param n_cores Number of cores.
#' @param fdr_thresh FDR threshold.
#' @param stability_thresh Frequency threshold to keep an edge (e.g., 0.8 = 80%).
#' @return List with edge table and network objects.
run_differential_network <- function(mat_ctrl, mat_case, n_boot = 100, n_perm = 1000, 
                                     seed = 123, n_cores = 1, fdr_thresh = 0.1, 
                                     stability_thresh = 0.8) {
  
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)
  requireNamespace("foreach", quietly = TRUE)
  requireNamespace("corpcor", quietly = TRUE)
  
  message(sprintf("      [DiffNet] Start: Boot=%d, Perm=%d, Cores=%d", n_boot, n_perm, n_cores))
  
  # Setup Cluster
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(corpcor) })
  # Esporta le funzioni necessarie ai worker
  parallel::clusterExport(cl, varlist = c("infer_network_pcor", "boot_worker_pcor"), envir = environment())
  parallel::clusterSetRNGStream(cl, seed)
  
  # --- STEP 1: BOOTSTRAP (Stability) ---
  message("      [DiffNet] 1/3 Running Bootstrap Stability...")
  
  # Helper interno per eseguire il bootstrap su una matrice
  run_boot_internal <- function(mat, n_b) {
    n_samp <- nrow(mat)
    res_list <- foreach::foreach(i = 1:n_b, .packages = "corpcor") %dopar% {
      boot_worker_pcor(mat, n_samp, lambda_val = NULL)
    }
    # Rimuovi NULL (iterazioni fallite)
    Filter(Negate(is.null), res_list)
  }
  
  boots_ctrl <- run_boot_internal(mat_ctrl, n_boot)
  boots_case <- run_boot_internal(mat_case, n_boot)
  
  # Compute edge stability
  agg_ctrl <- aggregate_boot_results(boots_ctrl, alpha = 0.05)
  agg_case <- aggregate_boot_results(boots_case, alpha = 0.05)
  
  check_boot_yield(agg_ctrl, n_boot, "Control")
  check_boot_yield(agg_case, n_boot, "Case")
  
  # Stability mask: Edge is counted for test correction if its stable in at least one group
  stable_mask <- (agg_ctrl$adj == 1 | agg_case$adj == 1)
  
  n_stable_edges <- sum(stable_mask[upper.tri(stable_mask)])
  message(sprintf("      [DiffNet] Found %d stable edges (union of groups).", n_stable_edges))
  
  if (n_stable_edges == 0) {
    parallel::stopCluster(cl)
    return(list(edges_table = data.frame(), networks = list()))
  }
  
  # --- STEP 2: OBSERVED DIFFERENCE ---
  obs_ctrl <- infer_network_pcor(mat_ctrl, fixed_lambda = NULL)
  obs_case <- infer_network_pcor(mat_case, fixed_lambda = NULL)
  obs_diff <- abs(obs_ctrl - obs_case)
  
  # --- STEP 3: PERMUTATION TEST (Only on Stable Edges logic) ---
  message("      [DiffNet] 2/3 Running Permutation Test...")
  
  pool_mat <- rbind(mat_ctrl, mat_case)
  n1 <- nrow(mat_ctrl)
  n2 <- nrow(mat_case)
  n_total <- n1 + n2
  
  null_diffs <- foreach::foreach(i = 1:n_perm, .packages = "corpcor") %dopar% {
    shuffled_idx <- sample(1:n_total)
    p1 <- pool_mat[shuffled_idx[1:n1], ]
    p2 <- pool_mat[shuffled_idx[(n1 + 1):n_total], ]
    
    r1 <- infer_network_pcor(p1, fixed_lambda = NULL)
    r2 <- infer_network_pcor(p2, fixed_lambda = NULL)
    abs(r1 - r2)
  }
  
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  # --- STEP 4: P-VALUE & FDR ---
  message("      [DiffNet] 3/3 calculating statistics...")
  
  nodes <- colnames(mat_ctrl)
  edges_idx <- which(upper.tri(stable_mask) & stable_mask, arr.ind = TRUE)
  
  results_list <- list()
  denom <- n_perm + 1
  
  if(nrow(edges_idx) > 0) {
    for(k in 1:nrow(edges_idx)) {
      i <- edges_idx[k,1]
      j <- edges_idx[k,2]
      
      obs_val <- obs_diff[i,j]
      null_vals <- sapply(null_diffs, function(m) m[i, j])
      
      p_val <- (sum(null_vals >= obs_val) + 1) / denom
      
      results_list[[k]] <- data.frame(
        Node1 = nodes[i],
        Node2 = nodes[j],
        Weight_Ctrl = obs_ctrl[i,j],
        Weight_Case = obs_case[i,j],
        Is_Stable_Ctrl = agg_ctrl$adj[i,j] == 1,
        Is_Stable_Case = agg_case$adj[i,j] == 1,
        Diff_Score = obs_val,
        P_Value = p_val,
        stringsAsFactors = FALSE
      )
    }
    
    res_df <- do.call(rbind, results_list)
    
    # FDR Correction (Only on stable edges)
    res_df$FDR <- stats::p.adjust(res_df$P_Value, method = "BH")
    res_df$Significant <- res_df$FDR < fdr_thresh
    
    res_df <- res_df[order(res_df$P_Value), ]
  } else {
    res_df <- data.frame()
  }
  
  return(list(
    edges_table = res_df,
    networks = list(ctrl = obs_ctrl, case = obs_case),
    stability = list(ctrl = agg_ctrl$adj, case = agg_case$adj)
  ))
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
