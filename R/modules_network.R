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
#' @param label_ctrl Label for control group (used for dynamic column naming).
#' @param label_case Label for case group (used for dynamic column naming).
#' @return List with edge table and network objects.
run_differential_network <- function(mat_ctrl, mat_case, n_boot = 100, n_perm = 1000, 
                                     seed = 123, n_cores = 1, fdr_thresh = 0.1, 
                                     stability_thresh = 0.8,
                                     label_ctrl = "Ctrl", label_case = "Case") {
  
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)
  requireNamespace("foreach", quietly = TRUE)
  requireNamespace("corpcor", quietly = TRUE)
  
  message(sprintf("      [DiffNet] Start: Boot=%d, Perm=%d, Cores=%d", n_boot, n_perm, n_cores))
  
  # Setup Cluster
  cl <- parallel::makeCluster(n_cores)
  # SAFETY: Ensure cluster is stopped even if code crashes
  on.exit({
    if(!is.null(cl)) {
      try(parallel::stopCluster(cl), silent=TRUE)
    }
  }, add = TRUE)
  
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(corpcor) })
  parallel::clusterExport(cl, varlist = c("infer_network_pcor", "boot_worker_pcor"), envir = environment())
  parallel::clusterSetRNGStream(cl, seed)
  
  # --- STEP 1: BOOTSTRAP (Stability) ---
  message("      [DiffNet] 1/3 Running Bootstrap Stability...")
  
  # Helper internal function
  run_boot_internal <- function(mat, n_b) {
    n_samp <- nrow(mat)
    res_list <- foreach::foreach(i = 1:n_b, .packages = "corpcor") %dopar% {
      boot_worker_pcor(mat, n_samp, lambda_val = NULL)
    }
    Filter(Negate(is.null), res_list)
  }
  
  boots_ctrl <- run_boot_internal(mat_ctrl, n_boot)
  boots_case <- run_boot_internal(mat_case, n_boot)
  
  # Compute edge stability
  agg_ctrl <- aggregate_boot_results(boots_ctrl, alpha = 0.05)
  agg_case <- aggregate_boot_results(boots_case, alpha = 0.05)
  
  check_boot_yield(agg_ctrl, n_boot, "Ctrl")
  check_boot_yield(agg_case, n_boot, "Case")
  
  # Stability mask
  stable_mask <- (agg_ctrl$adj == 1 | agg_case$adj == 1)
  
  n_stable_edges <- sum(stable_mask[upper.tri(stable_mask)])
  message(sprintf("      [DiffNet] Found %d stable edges (union of groups).", n_stable_edges))
  
  if (n_stable_edges == 0) {
    return(list(edges_table = data.frame(), networks = list()))
  }
  
  # --- STEP 2: OBSERVED DIFFERENCE ---
  obs_ctrl <- infer_network_pcor(mat_ctrl, fixed_lambda = NULL)
  obs_case <- infer_network_pcor(mat_case, fixed_lambda = NULL)
  obs_diff <- abs(obs_ctrl - obs_case)
  
  # --- STEP 3: PERMUTATION TEST ---
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
    
    res_df$FDR <- stats::p.adjust(res_df$P_Value, method = "BH")
    res_df$Significant <- res_df$FDR < fdr_thresh
    res_df <- res_df[order(res_df$P_Value), ]
    
    # --- DYNAMIC RENAMING ---
    colnames(res_df)[colnames(res_df) == "Weight_Ctrl"] <- paste0("Weight_", label_ctrl)
    colnames(res_df)[colnames(res_df) == "Weight_Case"] <- paste0("Weight_", label_case)
    colnames(res_df)[colnames(res_df) == "Is_Stable_Ctrl"] <- paste0("Is_Stable_", label_ctrl)
    colnames(res_df)[colnames(res_df) == "Is_Stable_Case"] <- paste0("Is_Stable_", label_case)
    
  } else {
    res_df <- data.frame()
  }
  
  return(list(
    edges_table = res_df,
    networks = list(ctrl = obs_ctrl, case = obs_case),
    stability = list(ctrl = agg_ctrl$adj, case = agg_case$adj)
  ))
}

#' @title Calculate Node Topology Metrics
#' @description 
#' Calculates centralities and topology metrics for a given adjacency matrix.
#' Returns a clean node-level dataframe.
#' 
#' @param adj_mat Adjacency matrix (0/1).
#' @return Dataframe with node metrics.
calculate_node_topology <- function(adj_mat) {
  
  requireNamespace("igraph", quietly = TRUE)
  
  if (sum(adj_mat) == 0) return(NULL)
  
  # Build Graph
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Calculate Metrics
  metrics <- data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    Cluster_Coeff = igraph::transitivity(g, type = "local"),
    stringsAsFactors = FALSE
  )
  
  # Order by Degree (Hubs first)
  metrics <- metrics[order(metrics$Degree, decreasing = TRUE), ]
  
  return(metrics)
}

#' @title Calculate Node Jaccard Rewiring
#' @description 
#' Quantifies how much a node's neighborhood changes between two networks.
#' Jaccard = Intersection(Neighbors) / Union(Neighbors).
#' 
#' @param adj_ctrl Adjacency matrix (0/1) for Control.
#' @param adj_case Adjacency matrix (0/1) for Case.
#' @return Dataframe with Jaccard index per node.
calculate_jaccard_rewiring <- function(adj_ctrl, adj_case) {
  
  # Ensure node alignment
  nodes <- intersect(rownames(adj_ctrl), rownames(adj_case))
  if (length(nodes) == 0) return(NULL)
  
  res_df <- data.frame(
    Node = nodes,
    Degree_Ctrl = 0,
    Degree_Case = 0,
    Jaccard_Index = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(nodes)) {
    n <- nodes[i]
    
    # Get neighbors (names of connected nodes)
    # Using which() != 0 is safer for float weights, though adj should be binary
    neigh_ctrl <- names(which(adj_ctrl[n, ] != 0))
    neigh_case <- names(which(adj_case[n, ] != 0))
    
    # Sets
    intersect_set <- intersect(neigh_ctrl, neigh_case)
    union_set <- union(neigh_ctrl, neigh_case)
    
    # Store Degrees
    res_df$Degree_Ctrl[i] <- length(neigh_ctrl)
    res_df$Degree_Case[i] <- length(neigh_case)
    
    # Calculate Jaccard
    if (length(union_set) > 0) {
      res_df$Jaccard_Index[i] <- length(intersect_set) / length(union_set)
    } else {
      # Both empty: Stable isolation. 
      # Arguably 1 (perfectly preserved state) or NA. 
      # We set to 1 to indicate "No rewiring occurred".
      res_df$Jaccard_Index[i] <- 1.0 
    }
  }
  
  # Add interpretation
  res_df$Rewiring_Status <- cut(res_df$Jaccard_Index, 
                                breaks = c(-Inf, 0.3, 0.7, Inf), 
                                labels = c("High_Rewiring", "Moderate", "Stable"))
  
  return(res_df)
}

#' @title Integrate Drivers and Topology (Hub-Driver Analysis)
#' @description 
#' Merges PLS-DA drivers with Network Topology metrics.
#' 
#' @param drivers_df Dataframe output from extract_plsda_loadings.
#' @param topo_df Dataframe output from calculate_node_topology (usually Case network).
#' @return Merged dataframe.
integrate_hub_drivers <- function(drivers_df, topo_df) {
  
  if (is.null(drivers_df) || nrow(drivers_df) == 0) return(NULL)
  if (is.null(topo_df) || nrow(topo_df) == 0) return(NULL)
  
  # Prepare Drivers: We need a single "Importance" metric.
  # If 'Importance' column exists (calculated in modules_multivariate), use it.
  # Otherwise, take max absolute weight.
  
  df_drv <- drivers_df
  if (!"Importance" %in% names(df_drv)) {
    # Fallback: Calculate max weight across available components
    w_cols <- grep("Weight|Comp", names(df_drv), value = TRUE)
    if (length(w_cols) > 0) {
      df_drv$Importance <- apply(df_drv[, w_cols, drop=FALSE], 1, function(x) max(abs(x), na.rm=TRUE))
    } else {
      return(NULL) # Cannot determine importance
    }
  }
  
  # Merge
  merged <- merge(df_drv, topo_df, by.x = "Marker", by.y = "Node")
  
  # Calculate Quartiles for Categorization
  imp_med <- median(merged$Importance, na.rm = TRUE)
  deg_med <- median(merged$Degree, na.rm = TRUE)
  
  merged$Role <- case_when(
    merged$Importance >= imp_med & merged$Degree >= deg_med ~ "Master_Regulator",
    merged$Importance >= imp_med & merged$Degree < deg_med  ~ "Solo_Driver",
    merged$Importance < imp_med & merged$Degree >= deg_med  ~ "Structural_Connector",
    TRUE ~ "Background"
  )
  
  return(merged)
}

#' @title Calculate Topology and Export Rich Cytoscape Table
#' @description 
#' Calculates node topology metrics and merges them into the edge list.
#' Creates a single "Rich Table" for Cytoscape (One-click import).
#' 
#' @param adj_mat Adjacency matrix (0/1).
#' @param weight_mat Partial correlation matrix.
#' @param group_label Label for the group (e.g., "Healthy", "LS").
#' @return A dataframe compatible with Cytoscape import.
get_rich_network_table <- function(adj_mat, weight_mat, group_label = "Unknown") {
  
  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # 1. Check if network is empty
  if (sum(adj_mat) == 0) return(NULL)
  
  # 2. Build Graph & Calculate Topology
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Calculate Metrics
  node_metrics <- data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )
  
  # 3. Build Edge List
  edge_list <- igraph::as_data_frame(g, what = "edges")
  
  # 4. Merge Logic (The "Rich" Table)
  # We join node metrics TWICE: once for Source, once for Target
  rich_table <- edge_list %>%
    rowwise() %>%
    mutate(
      Weight = weight_mat[from, to],
      Abs_Weight = abs(Weight),
      Sign = ifelse(Weight > 0, "Positive", "Negative"),
      Interaction = ifelse(Weight > 0, "Co-occurrence", "Mutual-Exclusion"),
      Group_ID = group_label
    ) %>%
    ungroup() %>%
    # Join Source Metrics
    left_join(node_metrics, by = c("from" = "Node")) %>%
    dplyr::rename(
      Source = from, 
      Source_Degree = Degree, 
      Source_Betweenness = Betweenness,
      Source_Closeness = Closeness
    ) %>%
    # Join Target Metrics
    left_join(node_metrics, by = c("to" = "Node")) %>%
    dplyr::rename(
      Target = to, 
      Target_Degree = Degree, 
      Target_Betweenness = Betweenness,
      Target_Closeness = Closeness
    ) %>%
    # Cleanup auxiliary columns from join if any
    select(-matches("Eigen_Centrality")) # Optional: keep or remove based on preference
  
  return(rich_table)
}
