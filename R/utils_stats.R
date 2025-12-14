# R/utils_stats.R
# ==============================================================================
# STATISTICAL INFERENCE ENGINE
# Description: Robust Partial Correlation, Bootstrap, Permutation
# ==============================================================================

library(corpcor)
library(doParallel)
library(foreach)

#' Convert Precision Matrix to Partial Correlation
#' @param prec_mat Inverse covariance matrix
prec2part <- function(prec_mat) {
  d <- diag(prec_mat)
  p <- ncol(prec_mat)
  part <- matrix(0, p, p)
  # Standardize to Correlation: -w_ij / sqrt(w_ii * w_jj)
  # Vectorized calculation for speed
  denom <- sqrt(outer(d, d))
  part <- -prec_mat / denom
  diag(part) <- 1
  return(part)
}

#' Single Estimate of Robust Partial Correlation
#' @param data Matrix (n x p)
#' @return Matrix (p x p)
estimate_pcor_robust <- function(data) {
  # Schafer-Strimmer Shrinkage
  fit <- corpcor::invcov.shrink(data, verbose = FALSE)
  return(prec2part(as.matrix(fit)))
}

#' Bootstrap Inference (Parallelized)
#' @param data Matrix (n x p)
#' @param B Integer, number of bootstraps
#' @param cl Parallel cluster object (optional)
#' @return List (adj_matrix, pcor_mean)
run_bootstrap_inference <- function(data, B = 1000, cl = NULL) {
  
  n <- nrow(data)
  p <- ncol(data)
  nodes <- colnames(data)
  
  # Handle cluster internally if not provided
  own_cluster <- FALSE
  if(is.null(cl)) {
    message("    [Stats] No cluster provided. Init local cluster...")
    n_cores <- parallel::detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    own_cluster <- TRUE
  }
  
  # Export dependencies to workers
  clusterEvalQ(cl, {
    library(corpcor)
  })
  
  # Parallel Loop
  boot_results <- foreach(b = 1:B, .combine = 'c', .packages='corpcor') %dopar% {
    idx <- sample(1:n, n, replace = TRUE)
    data_boot <- data[idx, ]
    
    tryCatch({
      fit <- corpcor::invcov.shrink(data_boot, verbose = FALSE)
      # Manual inversion to avoid function dependency issues on workers
      d <- diag(fit)
      p_mat <- -as.matrix(fit) / sqrt(outer(d, d))
      diag(p_mat) <- 1
      list(p_mat)
    }, error = function(e) return(NULL))
  }
  
  # Close cluster only if we created it
  if(own_cluster) stopCluster(cl)
  
  # Aggregation Logic
  valid_boots <- Filter(Negate(is.null), boot_results)
  n_valid <- length(valid_boots)
  
  message(sprintf("    [Stats] Valid Bootstraps: %d/%d", n_valid, B))
  
  # Compute CIs per edge
  adj_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  pcor_mean <- matrix(0, p, p, dimnames = list(nodes, nodes))
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      vals <- sapply(valid_boots, function(m) m[i,j])
      ci <- quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
      mu <- mean(vals, na.rm=TRUE)
      
      pcor_mean[i,j] <- mu; pcor_mean[j,i] <- mu
      
      # Significance: CI excludes 0
      if (ci[1] > 0 || ci[2] < 0) {
        adj_mat[i,j] <- 1; adj_mat[j,i] <- 1
      }
    }
  }
  
  return(list(adj = adj_mat, pcor_mean = pcor_mean))
}
