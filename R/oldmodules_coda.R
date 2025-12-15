# R/modules_coda.R
# ==============================================================================
# COMPOSITIONAL DATA ANALYSIS MODULE
# Description: Pure functions for Zero Replacement, Imputation, and CLR.
# Dependencies: zCompositions, robCompositions
# ==============================================================================

library(zCompositions)
library(robCompositions)
library(dplyr)

#' @title Zero Replacement via Count Zero Multiplicative (CZM) method
#' @description 
#' Replaces zeros in a compositional matrix using Bayesian-Multiplicative replacement.
#' Handles the "NA conflict" by temporarily masking NAs to allow cmultRepl to run.
#' 
#' @param mat A numeric matrix (samples x markers).
#' @return A matrix with zeros replaced by small estimated values, and NAs preserved.
replace_zeros_czm <- function(mat) {
  
  # Validation: Ensure input is a matrix
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  
  # Check if zeros exist
  if (!any(mat == 0, na.rm = TRUE)) {
    message("   [CoDa] No zeros detected. Skipping replacement.")
    return(mat)
  }
  
  message("   [CoDa] Zeros detected. Applying zCompositions::cmultRepl (CZM)...")
  
  # --- DEBUGGER FIX: HANDLE NAs BEFORE CMULTREPL ---
  # cmultRepl fails if it sees NAs. We perform a "Mask & Restore" strategy.
  
  # 1. Identify NA locations
  na_locs <- is.na(mat)
  has_na <- any(na_locs)
  
  # 2. Temporary Fill (if NAs exist)
  # We fill NAs with the column geometric mean (neutral value in CoDa) 
  # just to let cmultRepl focus on the Zeros without crashing.
  mat_temp <- mat
  
  if (has_na) {
    # Helper for geometric mean ignoring NA and <= 0
    gmean_safe <- function(x) {
      vals <- x[!is.na(x) & x > 0]
      if(length(vals) == 0) return(1e-6) # Fallback if col is empty/zeros
      return(exp(mean(log(vals))))
    }
    
    # Fill NAs column-wise
    for(j in 1:ncol(mat_temp)) {
      if(any(na_locs[,j])) {
        fill_val <- gmean_safe(mat_temp[,j])
        mat_temp[na_locs[,j], j] <- fill_val
      }
    }
  }
  
  # 3. Apply cmultRepl on the "complete" (masked) matrix
  tryCatch({
    clean_mat <- zCompositions::cmultRepl(
      X = mat_temp,
      label = 0,
      method = "CZM", 
      output = "p-counts",
      suppress.print = TRUE 
    )
    
    # 4. Restore NAs (The actual Imputation will be done by kNN later)
    if (has_na) {
      clean_mat[na_locs] <- NA
    }
    
    return(clean_mat)
    
  }, error = function(e) {
    stop(paste("Zero replacement failed:", e$message))
  })
}

#' @title Calculate Aitchison RMSD Error
#' @description Internal helper to calculate error between original and imputed values in log-ratio space.
calc_aitchison_error <- function(orig, imp, mask) {
  # Extract values where we artificially created NAs
  val_orig <- orig[mask]
  val_imp  <- imp[mask]
  
  # Avoid log(0) or log(neg) issues in error calc (safety net)
  val_imp[val_imp <= 0] <- min(val_orig[val_orig > 0])
  
  # Root Mean Square Difference of logs
  rmsd <- sqrt(mean((log(val_orig) - log(val_imp))^2))
  return(rmsd)
}

#' @title Optimize k for kNN Imputation (Cross-Validation)
#' @description 
#' Tests a range of k values by masking 5% of valid data and measuring imputation error.
#' 
#' @param mat A numeric matrix (no zeros).
#' @param k_seq A sequence of integers (e.g., 3:10).
#' @return A list containing: best_k (int), plot (ggplot object), results (dataframe).
optimize_knn_k <- function(mat, k_seq = 3:10) {
  
  message("[CoDa] Optimizing kNN parameter (k)...")
  
  # 1. Create a Test Set (Mask 5% of data)
  set.seed(123) # Reproducibility for the mask
  valid_indices <- which(!is.na(mat), arr.ind = TRUE)
  n_test <- floor(0.05 * nrow(valid_indices))
  
  if (n_test == 0) stop("Dataset too small for cross-validation.")
  
  test_subset_idx <- valid_indices[sample(1:nrow(valid_indices), n_test), ]
  
  # Create masked matrix
  mat_cv <- mat
  mat_cv[test_subset_idx] <- NA
  df_cv <- as.data.frame(mat_cv) # robCompositions needs DF
  
  errors <- numeric(length(k_seq))
  names(errors) <- k_seq
  
  # 2. Iterate over k
  for (i in seq_along(k_seq)) {
    k_curr <- k_seq[i]
    tryCatch({
      # Suppress output from impKNNa to keep console clean
      utils::capture.output({
        res <- robCompositions::impKNNa(df_cv, k = k_curr, method = "knn")
      })
      mat_imp <- as.matrix(res$xImp)
      
      # Calculate Error
      errors[i] <- calc_aitchison_error(mat, mat_imp, test_subset_idx)
    }, error = function(e) {
      errors[i] <- NA
      warning(sprintf("k=%d failed: %s", k_curr, e$message))
    })
  }
  
  # 3. Determine Winner
  results_df <- data.frame(k = k_seq, RMSE = errors)
  valid_res <- results_df %>% filter(!is.na(RMSE))
  
  if (nrow(valid_res) == 0) stop("All k optimizations failed.")
  
  best_k <- valid_res$k[which.min(valid_res$RMSE)]
  min_err <- min(valid_res$RMSE)
  
  message(sprintf("   -> Optimal k found: %d (RMSE: %.4f)", best_k, min_err))
  
  # 4. Generate Plot
  p <- ggplot(valid_res, aes(x = k, y = RMSE)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(size = 3, color = "steelblue") +
    geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
    labs(title = "kNN Imputation Optimization",
         subtitle = sprintf("Optimal k = %d (Minimizes Aitchison Error)", best_k),
         y = "Aitchison RMSE", x = "Neighbors (k)") +
    theme_bw()
  
  return(list(best_k = best_k, plot = p, results = results_df))
}

#' @title Robust Imputation for Missing Values (Wrapper)
#' @description Handles optimization logic internally if k is NULL.
impute_knn_aitchison <- function(mat, k = 3) {
  # Logic already defined in previous turn, just ensuring it accepts 'k'
  # This function calls robCompositions::impKNNa directly
  res_imp <- robCompositions::impKNNa(as.data.frame(mat), k = k, method = "knn")
  return(as.matrix(res_imp$xImp))
}

#' @title Centered Log-Ratio Transformation (CLR)
#' @description 
#' Projects the compositional data from the simplex to real space (Euclidean).
#' 
#' @param mat A numeric matrix (positive values only).
#' @return A CLR-transformed matrix.
transform_clr <- function(mat) {
  
  if (any(mat <= 0, na.rm = TRUE)) stop("CLR input must be strictly positive.")
  if (any(is.na(mat))) stop("CLR input must not contain NAs.")
  
  message("   [CoDa] Applying CLR transformation...")
  
  # Define Geometric Mean function
  gmean <- function(x) exp(mean(log(x)))
  
  # Apply CLR row-wise
  # formula: clr(x) = [ln(x1/g(x)), ..., ln(xp/g(x))]
  clr_mat <- t(apply(mat, 1, function(x) log(x / gmean(x))))
  
  return(clr_mat)
}